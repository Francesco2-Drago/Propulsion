function P = performancePoint(of, p_bar, thermo, engine)
    % performancePoint
    % Nozzle performance at one operating point (O/F, p_c). This is the single
    % thermochemistry-to-thrust evaluator: phase A, phase B and the acceptance
    % tests all call it, so a shape ranked in phase A and the same shape sized
    % in phase B are scored with identical physics.
    %
    % Effective exhaust velocity, with the throat area eliminated through
    % At = mdot*cstar/p_c:
    %   F = mdot*ue + (pe - pamb)*eps*At = mdot*[ue + (pe - pamb)*eps*cstar/p_c]
    % so c_eff = ue + (pe - pamb)*eps*cstar/p_c and F = mdot*c_eff, which is
    % INDEPENDENT of At at fixed p_c. With the assignment's pamb = 0 this is
    % c_eff = ue + eps*(pe/pc)*cstar.
    %
    % Inputs are clamped to the edges of the CEA lookup, as everywhere else in
    % the model: outside its domain the interpolant extrapolates.
    % INPUT
    %   of     / double / 1xM   Mixture ratio O/F [-] (vector allowed)
    %   p_bar  / double / 1x1   Chamber pressure [bar]
    %   thermo / struct / 1x1   Tc_fun, R_fun, k_fun, handles of (of, p_bar),
    %                           plus the lookup domain of_min, of_max,
    %                           p_bar_min, p_bar_max
    %   engine / struct / 1x1   eps [-] nozzle area ratio, pamb [Pa]
    % OUTPUT
    %   P      / struct / 1x1   Fields sized 1xM: Tc [K], R [J/kgK], gamma [-],
    %                           Me [-], pe_pc [-], pe [Pa], Te [K], ue [m/s],
    %                           cstar [m/s], c_eff [m/s], Isp [s]

    g0 = 9.80665;                              % [m/s2]

    of = of(:).';                              % row vector
    m = numel(of);

    % Clamp to the lookup domain
    of_c = min(max(of, thermo.of_min), thermo.of_max);
    p_c = min(max(p_bar, thermo.p_bar_min), thermo.p_bar_max);
    pc_Pa = p_c * 1e5;                         % [Pa]

    [Tc, R, gam, Me, pe_pc, pe, Te, ue, cstar, c_eff, Isp] = deal(nan(1, m));

    for i = 1:m
        % Chamber thermochemistry
        Tc(i) = thermo.Tc_fun(of_c(i), p_c);           % [K]
        R(i) = thermo.R_fun(of_c(i), p_c);             % [J/kgK]
        gam(i) = thermo.k_fun(of_c(i), p_c);           % [-]

        % Characteristic velocity, through the Vandenkerckhove function
        Gamma = sqrt(gam(i)) * (2/(gam(i) + 1))^((gam(i) + 1)/(2*(gam(i) - 1)));
        cstar(i) = sqrt(R(i)*Tc(i)) / Gamma;           % [m/s]

        % Supersonic branch of the area-Mach relation at the fixed area ratio
        Me(i) = supersonic_mach_from_area_ratio(engine.eps, gam(i));

        % Isentropic expansion to the exit section
        stag = 1 + 0.5*(gam(i) - 1)*Me(i)^2;
        pe_pc(i) = stag^(-gam(i)/(gam(i) - 1));        % [-]
        pe(i) = pc_Pa * pe_pc(i);                      % [Pa]
        Te(i) = Tc(i) / stag;                          % [K]
        ue(i) = Me(i) * sqrt(gam(i)*R(i)*Te(i));       % [m/s]

        % Effective exhaust velocity and specific impulse
        c_eff(i) = ue(i) + (pe(i) - engine.pamb)*engine.eps*cstar(i)/pc_Pa;
        Isp(i) = c_eff(i) / g0;                        % [s]
    end

    P = struct("Tc", Tc, "R", R, "gamma", gam, "Me", Me, "pe_pc", pe_pc, ...
        "pe", pe, "Te", Te, "ue", ue, "cstar", cstar, "c_eff", c_eff, "Isp", Isp);
end

%% FUNCTIONS

function Me = supersonic_mach_from_area_ratio(area_ratio, k)
    % supersonic_mach_from_area_ratio
    % Exit Mach number on the supersonic branch of the area-Mach relation.
    % INPUT
    %   area_ratio / double / 1x1   Ae/At [-]
    %   k          / double / 1x1   Ratio of specific heats [-]
    % OUTPUT
    %   Me         / double / 1x1   Exit Mach number [-]

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
    % area_mach_ratio
    % Area ratio associated with a Mach number (isentropic, one-dimensional).
    % INPUT
    %   M       / double / 1x1   Mach number [-]
    %   k       / double / 1x1   Ratio of specific heats [-]
    % OUTPUT
    %   eps_val / double / 1x1   Area ratio A/At [-]

    eps_val = 1/M * sqrt(((1 + 0.5*(k - 1)*(M^2)) / (0.5*(k + 1)))^((k + 1)/(k - 1)));
end
