function M = shapeMerit(lookupN, Isp_of, cfg)
    % shapeMerit
    % Scale-free phase A evaluator: score a normalized grain shape (R_c = 1) by
    % the specific impulse it delivers ON THE LOADED MASS, each shape being
    % judged at its own best O/F level.
    %
    % Burnback is autonomous, db/dt = a*(mdot_ox/Ap)^n does not contain p_c, so
    % time is a quadrature and the whole O/F history factorizes as
    %   O/F(b~) = lambda * Phi~(b~),   Phi~ = Ap~^n/Pb~
    % where lambda = mdot_ox^(1-n) R_c^(2n-1)/(rho_f a L) carries every bit of
    % scale information. In phase A lambda is a free internal parameter: it is
    % swept, and the shape is scored at the level that suits it best. Nothing
    % here knows the size of the engine.
    %
    % THE OBJECTIVE IS COMPUTED DIRECTLY, NOT THROUGH SIGMA:
    %   R_c      = (t_b a mdot_ox^n / I~(b~_end))^(1/(2n+1))          from C6
    %   L        = mdot_ox^(1-n) R_c^(2n-1)/(rho_f a lambda)          from lambda
    %   m_load   = mdot_ox t_b + rho_f L (pi R_c^2 - Ap(0))
    %   Isp_load = I_tot/(g0 m_load)
    % m_load does NOT contain b~_end at all, and I_tot depends on it for less
    % than 0.2 % (the tail carries ~0.3 % of the mass at the lowest thrust), so
    % the objective is structurally insensitive to exactly where the burn is
    % declared over. Going through sigma instead made the answer hinge on that
    % choice: sigma moved from 0.0028 to 0.0303 depending on the termination
    % criterion, worth 3.4 s of objective, eleven times the spread of the whole
    % top ten. The decomposition
    %   Isp_load = Isp_med * (OF_med + 1)/(OF_med + 1/(1 - sigma))
    % stays exact and is still reported, but as a DIAGNOSTIC that splits the
    % loss into drift and sliver, never as the path of calculation.
    %
    % Time weighting: dt = Ap~^n db~/(a mdot_ox^n), so every time average over
    % the burn is an integral in b~ weighted by Ap~^n.
    %
    % C4 - THE OXIDIZER FLOW IS COMPUTED, NOT QUANTIFIED OVER A BOX.
    % mdot_ox is not free: C7 nails it down. From the thrust requirement,
    %   mdot_ox = F_target / mean_t[ (1 + 1/O_F(b~)) * c_eff(O_F(b~)) ]
    % with the same Ap~^n db~ time weight as every other average here. Every
    % term is available in phase A: O_F(b~) = lambda*Phi~(b~) from the sweep,
    % c_eff from the CEA interpolants at p_c = 20 bar. F_target is a delivery
    % requirement, not a phase B result, so this is not hardcoding.
    % With the flow known, G_ox(0) follows exactly from the 2.4 factorization
    % and C10 can be applied HERE, on the real flow, instead of asking whether
    % some flow in a search box might satisfy it. The existential form approved
    % shapes that only worked at mdot_ox = 30 kg/s, a flow the thrust
    % requirement forbids, and phase A parked on that corner.
    % END OF BURN: C5, C12 and burnout, whichever comes first. C12 needs R_c to
    % turn its 3 mm into a normalized length, and R_c is available here because
    % mdot_ox is, so the two are resolved together in two passes (I~ moves by
    % under 1 %, hence R_c by under 0.4 %). Phase A and phase B therefore stop
    % the burn at the same place, and the thermochemistry is clamped at of_max
    % rather than extrapolated past the table.
    % INPUT
    %   lookupN / struct / 1x1   Normalized geometry lookup at R_c = 1, fields
    %                            b [-], Ap [-], perim [-]
    %   Isp_of  / handle / 1x1   Isp [s] as a function of O/F [-], at the phase A
    %                            reference chamber pressure
    %   cfg     / struct / 1x1   n, a [-, m/s] regression law, rho_f [kg/m3],
    %                            of_min, of_max [-] (C5), F_target [N] (C7),
    %                            thrust_mode, t_b [s] (C6), Gox_lo, Gox_hi
    %                            [kg/m2 s] (C4), web_min [m] (C12), r_out [-]
    %                            outer radius ratio of the grain, ceff_of
    %                            [handle], g0, and optionally n_lambda (80)
    % OUTPUT
    %   M       / struct / 1x1   ok, merit = Isp_load [s], lambda [-], Isp_med
    %                            [s], OF_med [-], sigma [-], drift, drift_full
    %                            [-], Gtilde [-], mdot_ox [kg/s], Gox0
    %                            [kg/m2 s], R_c [m], L [m], m_load [kg], I_tot
    %                            [Ns], b_end [-], Phi0 [-], of0, of_end [-]

    M = struct("ok", false, "merit", -Inf, "lambda", NaN, "Isp_med", NaN, ...
        "OF_med", NaN, "sigma", NaN, "drift", NaN, "drift_full", NaN, ...
        "Gtilde", NaN, "mdot_ox", NaN, "Gox0", NaN, "R_c", NaN, "L", NaN, ...
        "m_load", NaN, "I_tot", NaN, ...
        "b_end", NaN, "Phi0", NaN, "of0", NaN, "of_end", NaN);

    n = cfg.n;
    if ~isfield(cfg, "n_lambda") || isempty(cfg.n_lambda)
        cfg.n_lambda = 80;
    end

    b = lookupN.b(:);
    Ap = lookupN.Ap(:);
    Pb = lookupN.perim(:);
    if numel(b) < 3 || any(~isfinite(Ap)) || any(~isfinite(Pb))
        return
    end

    % Shape function and its initial value
    Phi = Ap.^n ./ max(Pb, realmin);
    if Phi(1) <= 0 || ~isfinite(Phi(1))
        return
    end
    M.Phi0 = Phi(1);

    % Sweep the O/F level. The bounds put the INITIAL O/F exactly on the two
    % edges of the CEA domain (C5), so no admissible level is left out.
    lambda_lo = cfg.of_min / Phi(1);
    lambda_hi = cfg.of_max / Phi(1);
    lambda_grid = logspace(log10(lambda_lo), log10(lambda_hi), cfg.n_lambda);

    merit_grid = -inf(size(lambda_grid));
    for i = 1:numel(lambda_grid)
        E = evaluate_level(lambda_grid(i), b, Ap, Pb, Phi, Isp_of, cfg);
        if E.ok
            merit_grid(i) = E.Isp_load;
        end
    end

    [best, ibest] = max(merit_grid);
    if ~isfinite(best)
        return
    end

    % Refine inside the bracket around the best grid point: the merit has to be
    % smooth in the shape variables or the direct search chases grid noise
    lo = lambda_grid(max(1, ibest - 1));
    hi = lambda_grid(min(numel(lambda_grid), ibest + 1));
    if hi > lo
        neg = @(lam) -level_merit(lam, b, Ap, Pb, Phi, Isp_of, cfg);
        lam_opt = fminbnd(neg, lo, hi, optimset("Display", "off", "TolX", 1e-4));
    else
        lam_opt = lambda_grid(ibest);
    end

    % Report at the winning level, falling back to the grid point if the
    % refinement landed somewhere worse
    E = evaluate_level(lam_opt, b, Ap, Pb, Phi, Isp_of, cfg);
    if ~E.ok || E.Isp_load < best
        E = evaluate_level(lambda_grid(ibest), b, Ap, Pb, Phi, Isp_of, cfg);
        lam_opt = lambda_grid(ibest);
    end
    if ~E.ok
        return
    end

    M.ok = true;
    M.merit = E.Isp_load;
    M.lambda = lam_opt;
    M.Isp_med = E.Isp_med;
    M.OF_med = E.OF_med;
    M.sigma = E.sigma;
    M.drift = E.drift;
    M.drift_full = E.drift_full;
    M.Gtilde = E.Gtilde;
    M.mdot_ox = E.mdot_ox;
    M.Gox0 = E.Gox0;
    M.R_c = E.R_c;
    M.L = E.L;
    M.m_load = E.m_load;
    M.I_tot = E.I_tot;
    M.b_end = E.b_end;
    M.of0 = E.of0;
    M.of_end = E.of_end;
end

%% FUNCTIONS

function v = level_merit(lambda, b, Ap, Pb, Phi, Isp_of, cfg)
    % level_merit
    % Scalar wrapper of evaluate_level for the one-dimensional refinement.
    % INPUT
    %   lambda / double / 1x1   O/F level [-]
    %   b, Ap, Pb, Phi / double / Kx1   Normalized burnback lookup and Phi~
    %   Isp_of / handle / 1x1   Isp [s] as a function of O/F [-]
    %   cfg    / struct / 1x1   Evaluator configuration
    % OUTPUT
    %   v      / double / 1x1   Isp on the loaded mass [s], -Inf if invalid

    E = evaluate_level(lambda, b, Ap, Pb, Phi, Isp_of, cfg);
    if E.ok
        v = E.Isp_load;
    else
        v = -Inf;
    end
end

function E = evaluate_level(lambda, b, Ap, Pb, Phi, Isp_of, cfg)
    % evaluate_level
    % Evaluate one O/F level: find where the burn ends, size the engine from the
    % thrust requirement, and take the objective DIRECTLY from the total impulse
    % over the loaded mass.
    % INPUT
    %   lambda / double / 1x1   O/F level, O/F(b~) = lambda*Phi~(b~) [-]
    %   b, Ap, Pb, Phi / double / Kx1   Normalized burnback lookup and Phi~
    %   Isp_of / handle / 1x1   Isp [s] as a function of O/F [-]
    %   cfg    / struct / 1x1   See shapeMerit
    % OUTPUT
    %   E      / struct / 1x1   ok, Isp_load, Isp_med, OF_med, sigma, drift,
    %                           drift_full, Gtilde, mdot_ox, Gox0, R_c, L,
    %                           m_load, I_tot, b_end, of0, of_end

    E = struct("ok", false, "Isp_load", -Inf, "Isp_med", NaN, "OF_med", NaN, ...
        "sigma", NaN, "drift", NaN, "drift_full", NaN, "Gtilde", NaN, ...
        "mdot_ox", NaN, "Gox0", NaN, "R_c", NaN, "L", NaN, "m_load", NaN, ...
        "I_tot", NaN, "b_end", NaN, "of0", NaN, "of_end", NaN);

    n = cfg.n;
    of_full = lambda * Phi;
    if of_full(1) < cfg.of_min
        return                       % C5 violated from the first instant
    end

    % ------------ END OF BURN AND SIZE, RESOLVED TOGETHER ------------
    % C12 caps the burnback at the web left over the casing, (1 - r_out) minus
    % the 3 mm liner, but 3 mm is a physical length and this is a normalized
    % shape: it takes R_c, which takes mdot_ox, which takes the truncation. Two
    % passes close the loop, since I~ moves by under 1 % between them and R_c
    % goes as I~^(-1/(2n+1)), i.e. under 0.4 %.
    b_stop = Inf;
    for pass = 1:2
        [bc, Apc, Pbc, ofc] = truncate_burn(b, Ap, Pb, of_full, cfg.of_max, b_stop);
        if numel(bc) < 3
            return
        end

        w_t = Apc.^n;                                  % dt is proportional to this
        I_burn = trapz(bc, w_t);                       % I~(b~_end)
        if ~(I_burn > 0)
            return
        end

        % Thrust requirement C7 -> oxidizer flow. Clamped to the CEA domain
        % rather than extrapolated: outside it the interpolant means nothing.
        ofc_th = min(max(ofc, cfg.of_min), cfg.of_max);
        ceff = cfg.ceff_of(ofc_th);
        if any(~isfinite(ceff)) || any(ceff <= 0)
            return
        end
        if cfg.thrust_mode == "initial"
            thrust_per_flow = (1 + 1/ofc_th(1)) * ceff(1);
        else
            thrust_per_flow = trapz(bc, (1 + 1 ./ ofc_th) .* ceff .* w_t) / I_burn;
        end
        if ~(thrust_per_flow > 0)
            return
        end
        mdot_ox = cfg.F_target / thrust_per_flow;              % [kg/s]

        % C6 -> casing radius
        R_c = (cfg.t_b * cfg.a * mdot_ox^n / I_burn)^(1/(2*n + 1));   % [m]

        b_stop_new = (1 - cfg.r_out) - cfg.web_min/R_c;        % [-] C12
        if abs(b_stop_new - b_stop) < 1e-9
            break
        end
        b_stop = b_stop_new;
    end
    if ~(b_stop > 0)
        return                       % the web is thinner than the liner
    end

    % ------------ SIZE AND OBJECTIVE, DIRECTLY ------------
    % lambda = mdot_ox^(1-n) R_c^(2n-1)/(rho_f a L)  ->  L
    L = mdot_ox^(1 - n) * R_c^(2*n - 1) / (cfg.rho_f * cfg.a * lambda);   % [m]

    % Loaded mass: oxidizer plus ALL the fuel put on board, sliver included.
    % Note this does not contain b~_end: it is the whole annulus between the
    % initial port and the casing.
    m_load = mdot_ox*cfg.t_b + cfg.rho_f * L * R_c^2 * (pi - Ap(1));      % [kg]
    if ~(m_load > 0)
        return
    end

    % Total impulse. dt = R_c^(2n+1) Ap~^n db~/(a mdot_ox^n) and
    % F = mdot_ox (1 + 1/O_F) c_eff, so
    I_tot = mdot_ox^(1 - n) * R_c^(2*n + 1) / cfg.a * ...
        trapz(bc, (1 + 1 ./ ofc_th) .* ceff .* w_t);                     % [Ns]

    E.Isp_load = I_tot / (cfg.g0 * m_load);                              % [s]

    % ------------ C4: THE INITIAL FLUX, ON THE REAL FLOW ------------
    E.Gtilde = I_burn^(2/(2*n + 1)) / Ap(1);
    E.Gox0 = mdot_ox^(1/(2*n + 1)) / (cfg.t_b*cfg.a)^(2/(2*n + 1)) * E.Gtilde;
    if E.Gox0 < cfg.Gox_lo || E.Gox0 > cfg.Gox_hi
        return
    end

    % ------------ DIAGNOSTICS ------------
    % The decomposition Isp_load = Isp_med*(OF_med+1)/(OF_med+1/(1-sigma)) is
    % exact and is checked by the acceptance tests, but it is NOT how Isp_load
    % was obtained above.
    den_f = trapz(bc, Apc.^n ./ max(Phi_of(Apc, Pbc, n), realmin));
    if den_f > 0
        E.OF_med = lambda * I_burn / den_f;
    end

    w_m = (1 + 1 ./ ofc_th) .* w_t;
    Isp_vals = Isp_of(ofc_th);
    if all(isfinite(Isp_vals))
        E.Isp_med = trapz(bc, Isp_vals .* w_m) / trapz(bc, w_m);
    end

    fuel_loaded = pi - Ap(1);                     % casing area at R_c = 1 [-]
    if fuel_loaded > 0
        E.sigma = min(max((pi - Apc(end)) / fuel_loaded, 0), 1 - 1e-9);
    end

    % TWO drifts. The first is over the part of the burn that actually carries
    % mass and is the one to design against. The second runs to burnout: for a
    % star it explodes, because the tips reach the casing before the flats and
    % the burning perimeter collapses over the last ~1.4 % of the web. That is
    % not noise, it is the reason a star is a poor fit for this mission, and it
    % belongs in the discussion of point (ii).
    E.drift = max(ofc) / max(min(ofc), realmin);
    E.drift_full = max(of_full) / max(min(of_full), realmin);

    E.mdot_ox = mdot_ox;
    E.R_c = R_c;
    E.L = L;
    E.m_load = m_load;
    E.I_tot = I_tot;
    E.b_end = bc(end);
    E.of0 = ofc(1);
    E.of_end = ofc(end);
    E.ok = isfinite(E.Isp_load);
end

function Phi = Phi_of(Ap, Pb, n)
    % Phi_of
    % Shape function Phi~ = Ap~^n/Pb~ on an arbitrary sample set.
    % INPUT
    %   Ap / double / Kx1   Normalized port area [-]
    %   Pb / double / Kx1   Normalized burning perimeter [-]
    %   n  / double / 1x1   Regression exponent [-]
    % OUTPUT
    %   Phi / double / Kx1  Shape function [-]

    Phi = Ap.^n ./ max(Pb, realmin);
end

function [bc, Apc, Pbc, ofc] = truncate_burn(b, Ap, Pb, of, of_max, b_stop)
    % truncate_burn
    % Cut the burnback history where the burn actually ends: at the first
    % crossing of the upper CEA bound (C5), at the residual-web stop (C12), or
    % at burnout, whichever comes first. Both cuts are interpolated rather than
    % snapped to a lookup level, because a b~_end that jumps by one sample makes
    % the merit discontinuous in the shape variables, and the direct search
    % would happily optimize the jump.
    % INPUT
    %   b, Ap, Pb, of / double / Kx1   Burnback, port area, perimeter, O/F
    %   of_max        / double / 1x1   Upper O/F bound (C5) [-]
    %   b_stop        / double / 1x1   Residual-web stop (C12) [-], Inf when not
    %                                  yet known (first pass)
    % OUTPUT
    %   bc, Apc, Pbc, ofc / double / Jx1   Truncated history

    % C5: interpolated crossing of the CEA ceiling
    k = find(of > of_max, 1, "first");
    if isempty(k)
        b_c5 = Inf;
    elseif k == 1
        bc = b(1); Apc = Ap(1); Pbc = Pb(1); ofc = of(1);
        return
    else
        f = (of_max - of(k-1)) / (of(k) - of(k-1));
        b_c5 = b(k-1) + f*(b(k) - b(k-1));
    end

    b_end = min([b_c5, b_stop, b(end)]);
    if ~(b_end > b(1))
        bc = b(1); Apc = Ap(1); Pbc = Pb(1); ofc = of(1);
        return
    end

    keep = b < b_end;
    bc = b(keep);
    Apc = Ap(keep);
    Pbc = Pb(keep);
    ofc = of(keep);

    % Close the history on the end point, interpolated when it falls between
    % two lookup levels
    if b_end < b(end)
        bc(end+1) = b_end;
        Apc(end+1) = interp1(b, Ap, b_end, "linear");
        Pbc(end+1) = interp1(b, Pb, b_end, "linear");
        ofc(end+1) = interp1(b, of, b_end, "linear");
    else
        bc(end+1) = b(end);
        Apc(end+1) = Ap(end);
        Pbc(end+1) = Pb(end);
        ofc(end+1) = of(end);
    end
end
