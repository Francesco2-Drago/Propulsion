function out = run_combustion(lookup, sim)
    % run_combustion
    % Exact combustion/burnback simulation (STEP 2-4 of the model) on a
    % PRECOMPUTED geometry lookup, i.e. WITHOUT rebuilding the MDF. Same physics
    % as run_combustion_case / combustionWrapped, so the sizing loop and the
    % final validation run use one consistent model. The burn ends where the
    % model's own events fire (O/F, pressure, burnout).
    % INPUT
    %   lookup / struct / 1x1   Geometry lookup with fields b, Ap, perim,
    %                           dAp_db, dperim_db, b_max
    %   sim    / struct / 1x1   Simulation config (grain_length, throat_area,
    %                           mdot_ox, fuel props, thermo handles + derivatives,
    %                           of/p bounds, time_output, fine_ode, eps, pamb)
    % OUTPUT
    %   out    / struct / 1x1   ok flag; on success: time, in_time and summary
    %                           metrics (burn_time, mean_thrust, mean_Isp, ...)

    g0 = 9.80665;
    out = struct("ok", false);

    try
        % Build the vars struct the physics helpers expect
        vars = struct();
        vars.geometry.grain_length = sim.grain_length;
        vars.geometry.throat_area = sim.throat_area;
        vars.geometry.port_area = interp_geometry(0, lookup.b, lookup.Ap);
        vars.geometry.burning_perimeter = interp_geometry(0, lookup.b, lookup.perim);
        vars.fuel.a_rf = sim.a_rf;
        vars.fuel.n_rf = sim.n_rf;
        vars.fuel.rho_f = sim.rho_f;
        vars.combustion.mdot_ox = sim.mdot_ox;
        vars.combustion.Tc_fun = sim.Tc_fun;
        vars.combustion.R_fun = sim.R_fun;
        vars.combustion.k_fun = sim.k_fun;
        vars.combustion.of_min = sim.of_min;
        vars.combustion.of_max = sim.of_max;
        vars.combustion.p_min = sim.p_min;
        vars.combustion.p_max = sim.p_max;

        % STEP 2: steady-state initial pressure (guard the bracket)
        % A missing sign change is NOT infeasibility: it means the steady-state
        % pressure lies outside [p_min, p_max]. Report why, so a bracketing
        % failure is never mistaken for a rejected design.
        zmin = Z_chamber_stst(sim.p_min, vars);
        zmax = Z_chamber_stst(sim.p_max, vars);
        if ~isfinite(zmin) || ~isfinite(zmax)
            out.err = sprintf(['STEP 2: Z_chamber_stst not finite at the bracket ends ' ...
                '(Z(p_min) = %g, Z(p_max) = %g).'], zmin, zmax);
            return
        end
        if sign(zmin) == sign(zmax)
            if zmin > 0
                where = "above p_max (raise p_max or lower mdot_ox / grain area)";
            else
                where = "below p_min (lower p_min or raise mdot_ox / grain area)";
            end
            out.err = sprintf(['STEP 2: no steady-state root in [%.3g, %.3g] bar ' ...
                '(Z(p_min) = %g, Z(p_max) = %g, same sign): p_c sits %s.'], ...
                sim.p_min*1e-5, sim.p_max*1e-5, zmin, zmax, where);
            return
        end
        p0 = fzero(@(pc) Z_chamber_stst(pc, vars), [sim.p_min, sim.p_max]);

        % STEP 3: integrate chamber pressure and burnback with the model events
        if sim.fine_ode
            vars.combustion.dRTdp_fun_OF_p = sim.dRTdp_fun;
            vars.combustion.dRTdOF_fun_OF_p = sim.dRTdOF_fun;
        end
        y0 = [p0; 0];
        ode_opts = odeset("RelTol", 1e-7, "AbsTol", [1e-3, 1e-9], ...
            "Events", @(t,y) burnback_events(t, y, lookup, vars));
        [time, Y] = ode113(@(t,y) ode_pressure_burnback(t, y, vars, lookup, sim.fine_ode), ...
            sim.time_output, y0, ode_opts);

        in_time.time = time(:);
        in_time.pch = Y(:,1);
        in_time.burnback = Y(:,2);

        % STEP 4: reconstruct performance histories
        nt = numel(time);
        [Ap_h, perim_h, of_h, gox_h, Tch_h, pe_h] = deal(nan(nt,1));
        [mdotin_h, thrust_h, Isp_h] = deal(nan(nt,1));
        for j = 1:nt
            b = in_time.burnback(j);
            Ap = interp_geometry(b, lookup.b, lookup.Ap);
            perim = interp_geometry(b, lookup.b, lookup.perim);
            vars.geometry.port_area = Ap;
            vars.geometry.burning_perimeter = perim;

            pch = in_time.pch(j);
            % O/F and GOX straight from the geometry (no validating call, so a
            % point sitting exactly at the O/F event does not raise an error)
            GOX = sim.mdot_ox / Ap;
            rf = sim.a_rf * GOX^sim.n_rf;
            mdot_f = sim.rho_f * perim * sim.grain_length * rf;
            O_F = sim.mdot_ox / mdot_f;

            [ofT, pbT] = clamp_thermo_inputs(O_F, pch*1e-5, vars);
            Tch = sim.Tc_fun(ofT, pbT);
            k = sim.k_fun(ofT, pbT);
            R = sim.R_fun(ofT, pbT);
            Me = supersonic_mach_from_area_ratio(sim.eps, k);
            Te = Tch / (1 + 0.5*(k - 1)*Me^2);
            pe = pch / ((1 + 0.5*(k - 1)*Me^2)^(k/(k - 1)));
            ue = Me * sqrt(k*R*Te);
            mdot_in = sim.mdot_ox / O_F + sim.mdot_ox;
            K2 = k*((2/(k + 1))^((k + 1)/(k - 1)));
            cstar = sqrt(R*Tch/K2);
            mdot_nozzle = pch*sim.throat_area/cstar;
            thrust = mdot_nozzle*ue + (pe - sim.pamb)*sim.eps*sim.throat_area;

            Ap_h(j)=Ap; perim_h(j)=perim; of_h(j)=O_F; gox_h(j)=GOX;
            Tch_h(j)=Tch; pe_h(j)=pe; mdotin_h(j)=mdot_in;
            thrust_h(j)=thrust; Isp_h(j)=thrust/(mdot_nozzle*g0);
        end

        in_time.Ap = Ap_h; in_time.perim = perim_h;
        in_time.O_F = of_h; in_time.GOX = gox_h; in_time.Tch = Tch_h;
        in_time.pe = pe_h; in_time.mdot_in = mdotin_h;
        in_time.thrust = thrust_h; in_time.Isp = Isp_h;

        % Summary metrics (time-averaged over the actual burn)
        out.ok = true;
        out.time = in_time.time;
        out.in_time = in_time;
        out.burn_time = time(end) - time(1);
        out.total_impulse = trapz(in_time.time, thrust_h);
        out.m_prop = trapz(in_time.time, mdotin_h);
        out.mean_thrust = out.total_impulse / max(out.burn_time, realmin);
        out.mean_Isp = out.total_impulse / (g0 * max(out.m_prop, realmin));
        out.initial_thrust = thrust_h(1);
        out.max_pressure = max(in_time.pch);
        out.of0 = of_h(1); out.of_end = of_h(end);
        out.Gox0 = gox_h(1); out.Gox_end = gox_h(end);
    catch ME
        out.ok = false;
        out.err = ME.message;
    end
end

%% Local functions

function value = interp_geometry(b, b_lookup, values)
    b = min(max(b, b_lookup(1)), b_lookup(end));
    value = interp1(b_lookup, values, b, "linear");
end

function dy = ode_pressure_burnback(~, y, vars, lookup, fine_ode_boolean)
    pch = y(1);
    b = y(2);

    Ap = interp_geometry(b, lookup.b, lookup.Ap);
    perim = interp_geometry(b, lookup.b, lookup.perim);
    dAp_db = interp_geometry(b, lookup.b, lookup.dAp_db);
    dperim_db = interp_geometry(b, lookup.b, lookup.dperim_db);

    L = vars.geometry.grain_length;
    At = vars.geometry.throat_area;
    a_rf = vars.fuel.a_rf;
    n_rf = vars.fuel.n_rf;
    rho_f = vars.fuel.rho_f;
    mdot_ox = vars.combustion.mdot_ox;
    Tc_fun = vars.combustion.Tc_fun;
    R_fun = vars.combustion.R_fun;
    k_fun = vars.combustion.k_fun;

    Gox = mdot_ox / Ap;
    rf = a_rf * (Gox^n_rf);
    Ab = perim * L;
    mdot_f = rho_f * Ab * rf;
    mdot_in = mdot_f + mdot_ox;

    O_F = mdot_ox / mdot_f;
    pc_bar = pch * 1e-5;
    [O_F_thermo, pc_bar_thermo] = clamp_thermo_inputs(O_F, pc_bar, vars);
    Tc = Tc_fun(O_F_thermo, pc_bar_thermo);
    R = R_fun(O_F_thermo, pc_bar_thermo);
    k = k_fun(O_F_thermo, pc_bar_thermo);
    K2 = k*((2/(k + 1))^((k + 1)/(k - 1)));
    cstar = sqrt(R*Tc/K2);
    mdot_out = pch*At/cstar;
    rho_g = pch/(R*Tc);

    geometry.grain_length = L;
    geometry.port_area = Ap;
    geometry.d_area_area = (dAp_db / Ap) * rf;

    if fine_ode_boolean
        dRTdp = vars.combustion.dRTdp_fun_OF_p(O_F_thermo, pc_bar_thermo);
        dRTdOF = vars.combustion.dRTdOF_fun_OF_p(O_F_thermo, pc_bar_thermo);

        fine_ode.deltapc = pc_bar_thermo * dRTdp / (R*Tc);
        fine_ode.deltaOF = O_F_thermo * dRTdOF / (R*Tc);
        fine_ode.n_rf = n_rf;
        fine_ode.d_perim_perim = (dperim_db / max(perim, realmin)) * rf;

        dp_p = ode_chamber_pressure(mdot_in, mdot_out, rho_g, geometry, fine_ode);
    else
        dp_p = ode_chamber_pressure(mdot_in, mdot_out, rho_g, geometry);
    end

    dy = [dp_p*pch; rf];
end

function [value, isterminal, direction] = burnback_events(~, y, lookup, vars)
    pch = y(1);
    b = y(2);
    perim = interp_geometry(b, lookup.b, lookup.perim);
    O_F = mixture_ratio_from_geometry(b, lookup, vars);

    value = [
        lookup.b_max - b
        perim
        pch - vars.combustion.p_min
        vars.combustion.p_max - pch
        O_F - vars.combustion.of_min
        vars.combustion.of_max - O_F
    ];
    isterminal = ones(size(value));
    direction = [-1; -1; -1; -1; 0; 0];
end

function O_F = mixture_ratio_from_geometry(b, lookup, vars)
    Ap = interp_geometry(b, lookup.b, lookup.Ap);
    perim = interp_geometry(b, lookup.b, lookup.perim);
    Gox = vars.combustion.mdot_ox / Ap;
    rf = vars.fuel.a_rf * (Gox^vars.fuel.n_rf);
    mdot_f = vars.fuel.rho_f * perim * vars.geometry.grain_length * rf;
    O_F = vars.combustion.mdot_ox / mdot_f;
end

function [O_F_eval, p_bar_eval] = clamp_thermo_inputs(O_F, p_bar, vars)
    O_F_eval = min(max(O_F, vars.combustion.of_min), vars.combustion.of_max);
    p_bar_eval = min(max(p_bar, vars.combustion.p_min*1e-5), vars.combustion.p_max*1e-5);
end

function Me = supersonic_mach_from_area_ratio(area_ratio, k)
    if area_ratio < 1
        error("Nozzle area ratio must be >= 1.");
    end
    lo = 1 + 1e-8;
    hi = 2;
    while area_mach_ratio(hi, k) < area_ratio
        hi = hi * 1.5;
        if hi > 100
            error("Could not bracket supersonic exit Mach number.");
        end
    end
    Me = fzero(@(M) area_mach_ratio(M, k) - area_ratio, [lo, hi]);
end

function eps_val = area_mach_ratio(M, k)
    eps_val = 1/M * sqrt(((1 + 0.5*(k - 1)*(M^2)) / (0.5*(k + 1)))^((k + 1)/(k - 1)));
end
