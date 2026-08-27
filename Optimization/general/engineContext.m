function ctx = engineContext(thermo, params, C)
    % engineContext
    % Build everything the two optimization phases share for ONE oxidizer: the
    % keyed constraint table, the lookup domain the thermochemistry gets
    % clamped to, the nozzle data, and the Isp(O/F) and c_eff(O/F) curves at the
    % reference chamber pressure. Assembled once per oxidizer and passed down,
    % so the cost function never rebuilds an interpolant per evaluation.
    % INPUT
    %   thermo / struct / 1x1   Raw interpolant handles from get_thermo
    %   params / struct / 1x1   Combustion parameters (from combustion_params())
    %   C      / struct / 1xM   Constraint table from optimizationConstraints
    % OUTPUT
    %   ctx    / struct / 1x1   C, K              constraint table and keyed view
    %                           box               search box of the shape variables
    %                           oxidizer          name of the oxidizer in use
    %                           thermo            handles + lookup domain
    %                           engine            eps [-], pamb [Pa]
    %                           Isp_of, ceff_of   handles of O/F [-]
    %                           p_ref             reference pressure [bar]
    %                           merit_cfg         config for shapeMerit
    %                           h_min             C2 tip floor at R_c = 1 [-]

    g0 = 9.80665;                            % [m/s2]

    [~, ~, box] = optimizationConstraints();
    ctx.C = C;
    ctx.K = constraintsById(C);
    ctx.box = box;
    K = ctx.K;

    ctx.oxidizer = thermo.name;

    % Attach the domain the interpolants are valid on: O/F from C5, pressure
    % from the extent of the tables
    ctx.thermo = thermo;
    ctx.thermo.of_min = K.C5.lo;                          % [-]   C5
    ctx.thermo.of_max = K.C5.hi;                          % [-]   C5
    ctx.thermo.p_bar_min = params.combustion.p_min*1e-5;  % [bar] lookup extent
    ctx.thermo.p_bar_max = params.combustion.p_max*1e-5;  % [bar] lookup extent

    ctx.engine.eps = params.engine.eps;      % [-]
    ctx.engine.pamb = params.engine.pamb;    % [Pa]

    % Isp(O/F) and c_eff(O/F) at the reference pressure. Phase A runs at the C8
    % ceiling: the design leans on it, and the exact p_c history is phase B's
    % job. Tabulated once and interpolated, since the exit Mach number needs a
    % root find. c_eff is what turns a mixture ratio into the oxidizer flow the
    % thrust requirement implies
    ctx.p_ref = K.C8.lo * 1e-5;              % [bar]  C8
    of_grid = linspace(K.C5.lo, K.C5.hi, 400);
    P = performancePoint(of_grid, ctx.p_ref, ctx.thermo, ctx.engine);

    Isp_interp = griddedInterpolant(of_grid, P.Isp, "linear", "nearest");
    ceff_interp = griddedInterpolant(of_grid, P.c_eff, "linear", "nearest");
    ctx.Isp_of = @(of) Isp_interp(of);
    ctx.ceff_of = @(of) ceff_interp(of);

    % O/F of peak Isp and the corresponding cstar: only used as starting
    % guesses for the phase B length and throat, never as a constraint
    [~, i_best] = max(P.Isp);
    ctx.of_best = of_grid(i_best);           % [-]
    ctx.cstar_best = P.cstar(i_best);        % [m/s]

    % Phase A merit configuration
    ctx.merit_cfg.n = params.fuel.n_rf;       % [-]
    ctx.merit_cfg.a = params.fuel.a_rf;       % [m/s]/[kg/(m2 s)]^n
    ctx.merit_cfg.rho_f = params.fuel.rho_f;  % [kg/m3]
    ctx.merit_cfg.of_min = K.C5.lo;           % [-]   C5
    ctx.merit_cfg.of_max = K.C5.hi;           % [-]   C5
    ctx.merit_cfg.n_lambda = 80;              % [-]   O/F levels swept
    ctx.merit_cfg.g0 = g0;                    % [m/s2]
    ctx.merit_cfg.web_min = K.C12.lo;         % [m]   C12, residual web
    ctx.merit_cfg.r_out = NaN;                % [-]   outer radius ratio of the
                                              %       grain, set per shape by
                                              %       the caller

    % C4 is C10 anticipated in phase A, evaluated on the flow the thrust
    % requirement implies. The margin is there because phase A works at a
    % constant p_c = 20 bar and so gets mdot_ox wrong by a few percent: without
    % it, shapes that scrape through in A fail in B.
    ctx.merit_cfg.F_target = K.C7.lo;             % [N]   C7
    ctx.merit_cfg.thrust_mode = K.C7.value;       % "mean" or "initial"
    ctx.merit_cfg.t_b = K.C6.lo;                  % [s]   C6
    ctx.merit_cfg.Gox_lo = K.C4.lo;               % [kg/m2 s]  C4 = C10 + margin
    ctx.merit_cfg.Gox_hi = K.C4.hi;               % [kg/m2 s]
    ctx.merit_cfg.ceff_of = ctx.ceff_of;
    ctx.merit_cfg.Isp_of = ctx.Isp_of;

    % C2 at R_c = 1: the tip must be at least K.C2.lo cells tall
    ctx.h_min = K.C2.lo / params.mdf.grid_divisions;   % [-]
end

function K = constraintsById(C)
    % constraintsById
    % Id-keyed view of the constraint table, so a module handed the struct
    % array can read K.C10.lo instead of a magic array index.
    % INPUT
    %   C / struct / 1xM   Constraint table from optimizationConstraints
    % OUTPUT
    %   K / struct / 1x1   One field per constraint id

    K = struct();
    for i = 1:numel(C)
        K.(char(C(i).id)) = C(i);
    end
end
