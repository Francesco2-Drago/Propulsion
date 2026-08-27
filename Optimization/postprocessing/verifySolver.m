function V = verifySolver(design, opts)
    % verifySolver
    % The three sweeps of the brief 7.1: take a design, turn ONE knob at a time
    % with everything else frozen, and plot the target each knob is supposed to
    % control. The point is not accuracy, it is MONOTONICITY: sizeEngine solves
    % C6, C7 and C8 by bisection on a bracket with a verified sign change, and
    % that is only safe if each target moves one way with its own knob.
    %
    %   t_burn   vs R_c   expected increasing, d ln t/d ln R_c ~ 2n+1 = 2.5
    %   F_mean   vs L     expected increasing (marginal on the fuel-rich branch)
    %   p_c,peak vs A_t   expected decreasing with slope exactly -1
    %
    % The last one is exact because mdot(b)*cstar(b) does not contain A_t at
    % fixed geometry, so p_c = max_b[mdot cstar]/A_t. Any departure from -1 is a
    % coupling through the thermochemistry, i.e. a measure of how strongly p_c
    % feeds back on the gas properties.
    %
    % The middle one is the delicate one: the brief warns that dF/dL changes
    % sign on the fuel-rich branch when d ln c_eff/d(O/F) exceeds s^2/(1+s) with
    % s = 1/(O/F). If it is not monotone over the range the solver explores, the
    % bracket search will fail exactly as it did before the step-size fix.
    % INPUT
    %   design / struct / 1x1   Output of phaseB, or a struct with fields S
    %                           (sized engine), shape and ctx. Omit to size a
    %                           default design here
    %   opts   / struct / 1x1   Optional: n_points (default 11), span (default
    %                           0.35, the half-width of each sweep as a
    %                           fraction), do_plots (default true), oxidizer
    %                           (default "O2(L)", only used when design is
    %                           omitted)
    % OUTPUT
    %   V      / struct / 1x1   One field per sweep (Rc, L, At), each with the
    %                           knob values, the target values and the fitted
    %                           log-log slope

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, "n_points"), opts.n_points = 11; end
    if ~isfield(opts, "span"), opts.span = 0.35; end
    if ~isfield(opts, "do_plots"), opts.do_plots = true; end
    if ~isfield(opts, "oxidizer"), opts.oxidizer = "O2(L)"; end

    repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(genpath(repoRoot));

    % ------------ THE DESIGN TO PROBE ------------
    if nargin < 1 || isempty(design)
        fprintf("No design given: sizing a default one on %s...\n", opts.oxidizer);
        design = default_design(opts.oxidizer, repoRoot);
    end
    S = design.S;
    ctx = design.ctx;
    lookupN = design.lookupN;
    cfg = design.cfg;

    fprintf("\n%s\n", repmat('=', 1, 72));
    fprintf(" SOLVER VERIFICATION - three sweeps around the design\n");
    fprintf("%s\n", repmat('-', 1, 72));
    fprintf(" R_c = %.4f m, L = %.4f m, A_t = %.5e m2, mdot_ox = %.3f kg/s\n", ...
        S.R_c, S.L, S.At, S.mdot_ox);
    fprintf(" targets: t_b = %.1f s, F = %.2f kN, p_c = %.2f bar\n", ...
        S.burn_time, S.mean_thrust*1e-3, S.p_peak*1e-5);
    fprintf("%s\n", repmat('=', 1, 72));

    n = cfg.n;
    f = linspace(1 - opts.span, 1 + opts.span, opts.n_points);

    % A burn that reaches the end of the integration horizon is CENSORED, not
    % measured: ode113 stops at tspan(end), so the reported time is the horizon
    % and not the physical burn time. Such points are shown but kept out of the
    % slope fit, or they would drag it towards zero.
    t_horizon = cfg.time_output(end);

    % ------------ SWEEP 1: burn time against casing radius (C6) ------------
    V.Rc = sweep("R_c [m]", "t_burn [s]", S.R_c*f, ...
        @(x) simulateEngine(x, S.L, S.At, S.mdot_ox, ctx, lookupN, cfg), ...
        @(r) r.burn_time, 2*n + 1, t_horizon);

    % ------------ SWEEP 2: mean thrust against grain length (C7) ------------
    V.L = sweep("L [m]", "F_mean [kN]", S.L*f, ...
        @(x) simulateEngine(S.R_c, x, S.At, S.mdot_ox, ctx, lookupN, cfg), ...
        @(r) r.mean_thrust*1e-3, NaN, Inf);

    % ------------ SWEEP 3: peak pressure against throat area (C8) ------------
    V.At = sweep("A_t [m2]", "p_c,peak [bar]", S.At*f, ...
        @(x) simulateEngine(S.R_c, S.L, x, S.mdot_ox, ctx, lookupN, cfg), ...
        @(r) r.max_pressure*1e-5, -1, Inf);

    % ------------ VERDICT ------------
    fprintf("\n%s\n", repmat('-', 1, 72));
    verdict(V.Rc, "t_burn vs R_c", +1, 2*n + 1, 0.15);
    verdict(V.L, "F_mean vs L", +1, NaN, NaN);
    verdict(V.At, "p_c,peak vs A_t", -1, -1, 0.05);
    fprintf("%s\n\n", repmat('=', 1, 72));

    if opts.do_plots
        plot_sweeps(V, S);
    end
end

%% FUNCTIONS

function s = sweep(knob_name, target_name, x_values, run_fun, read_fun, ...
        slope_expected, t_horizon)
    % sweep
    % Turn one knob over a range, read the corresponding target, and fit the
    % log-log slope on the points that were actually measured.
    % INPUT
    %   knob_name      / string / 1x1   Axis label of the knob
    %   target_name    / string / 1x1   Axis label of the target
    %   x_values       / double / 1xM   Knob values to try
    %   run_fun        / handle / 1x1   Knob value -> run_combustion output
    %   read_fun       / handle / 1x1   Run output -> target value
    %   slope_expected / double / 1x1   Expected d ln y/d ln x, NaN if unknown
    %   t_horizon      / double / 1x1   Integration horizon [s]; a burn that
    %                                   reaches it is censored. Inf to disable
    % OUTPUT
    %   s              / struct / 1x1   x, y, censored, slope, slope_expected,
    %                                   knob, target

    m = numel(x_values);
    y = nan(1, m);
    censored = false(1, m);
    fprintf("\n %-12s -> %-16s\n", knob_name, target_name);
    fprintf(" %14s %16s\n", knob_name, target_name);
    for i = 1:m
        r = run_fun(x_values(i));
        if isstruct(r) && isfield(r, "ok") && r.ok
            y(i) = read_fun(r);
            % The burn hit the end of tspan: the value is a lower bound
            censored(i) = isfinite(t_horizon) && ...
                r.burn_time >= t_horizon*(1 - 1e-6);
            if censored(i)
                fprintf(" %14.6g %16.4f   censored, burn reached the %g s horizon\n", ...
                    x_values(i), y(i), t_horizon);
            else
                fprintf(" %14.6g %16.4f\n", x_values(i), y(i));
            end
        else
            if isstruct(r) && isfield(r, "err")
                fprintf(" %14.6g %16s   %s\n", x_values(i), "-", string(r.err));
            else
                fprintf(" %14.6g %16s\n", x_values(i), "-");
            end
        end
    end

    good = isfinite(y) & y > 0 & x_values > 0 & ~censored;
    if nnz(good) >= 3
        p = polyfit(log(x_values(good)), log(y(good)), 1);
        slope = p(1);
    else
        slope = NaN;
    end

    s = struct("x", x_values, "y", y, "censored", censored, "slope", slope, ...
        "slope_expected", slope_expected, ...
        "knob", string(knob_name), "target", string(target_name));
end

function verdict(s, name, direction, slope_expected, slope_tol)
    % verdict
    % Report whether a sweep is monotone in the expected direction and, when
    % there is a predicted slope, whether it matches.
    % INPUT
    %   s              / struct / 1x1   Sweep record
    %   name           / string / 1x1   Human name
    %   direction      / double / 1x1   +1 increasing, -1 decreasing
    %   slope_expected / double / 1x1   Predicted log-log slope, NaN if none
    %   slope_tol      / double / 1x1   Relative tolerance on the slope
    % OUTPUT
    %   None (prints to stdout)

    % Monotonicity is judged on the measured points only: a censored value is a
    % lower bound and cannot contradict a trend
    good = isfinite(s.y) & ~s.censored;
    d = diff(s.y(good));
    if isempty(d)
        fprintf(" %-18s NO DATA\n", name);
        return
    end
    mono = all(direction*d > 0);
    if mono
        tag = "monotone";
    else
        n_wrong = nnz(direction*d <= 0);
        tag = sprintf("NOT MONOTONE (%d of %d steps go the wrong way)", ...
            n_wrong, numel(d));
    end

    if isfinite(slope_expected) && isfinite(s.slope)
        err = abs(s.slope/slope_expected - 1);
        if err <= slope_tol
            stag = sprintf("slope %.3f vs %.3f expected, ok", s.slope, slope_expected);
        else
            stag = sprintf("slope %.3f vs %.3f expected, OFF by %.1f%%", ...
                s.slope, slope_expected, 100*err);
        end
    else
        stag = sprintf("slope %.3f", s.slope);
    end

    fprintf(" %-18s %s, %s (%d/%d points evaluated)\n", ...
        name, tag, stag, nnz(good), numel(s.y));
end

function plot_sweeps(V, S)
    % plot_sweeps
    % The appendix figure of the brief 7.1: the three sweeps side by side, with
    % the design point marked.
    % INPUT
    %   V / struct / 1x1   Sweep results
    %   S / struct / 1x1   The sized engine
    % OUTPUT
    %   None (creates a figure)

    figure("Name", "Solver verification - the three sweeps", "Color", "w");
    tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");

    marks = {S.R_c, S.burn_time; S.L, S.mean_thrust*1e-3; S.At, S.p_peak*1e-5};
    names = ["Rc", "L", "At"];
    for i = 1:3
        s = V.(names(i));
        nexttile
        plot(s.x, s.y, "o-", "LineWidth", 1.5, "Color", [0 0.45 0.74]);
        hold on
        plot(marks{i,1}, marks{i,2}, "p", "MarkerSize", 13, ...
            "MarkerFaceColor", [0.85 0.33 0.10], "MarkerEdgeColor", "k");
        grid on
        xlabel(s.knob);
        ylabel(s.target);
        if isfinite(s.slope_expected)
            title(sprintf("slope %.2f (expected %.2f)", s.slope, s.slope_expected));
        else
            title(sprintf("slope %.2f", s.slope));
        end
    end
end

function design = default_design(OX, repoRoot)
    % default_design
    % Size a plain circular-port engine so verifySolver can run on its own,
    % without waiting for a full phaseA/phaseB.
    % INPUT
    %   OX       / string / 1x1   Oxidizer name
    %   repoRoot / string / 1x1   Repository root
    % OUTPUT
    %   design   / struct / 1x1   S, ctx, lookupN, cfg

    params = combustion_params();
    params.geometry.type = "cylinder";
    params.engine.ext_diameter = 2.0;         % R_c = 1
    params.mdf.grid_divisions = 350;
    params.mdf.perimeter_from_area = false;

    C = optimizationConstraints();
    K = constraintsById(C);
    thermo = get_thermo(OX, repoRoot);
    ctx = engineContext(thermo, params, C);

    % Phase A picks the port fraction and the oxidizer flow
    x1_grid = linspace(0.20, 0.40, 25);
    J = inf(size(x1_grid));
    for i = 1:numel(x1_grid)
        J(i) = cylinderCostFunction(x1_grid(i), params, ctx);
    end
    [~, ibest] = min(J);
    [~, info] = cylinderCostFunction(x1_grid(ibest), params, ctx);
    if ~info.feasible
        error("verifySolver:noDesign", ...
            "Could not find an admissible default design for %s.", OX);
    end

    lookupN = build_shape_lookup(struct("diameter", 2*x1_grid(ibest)), params);
    cfg = struct("n", params.fuel.n_rf, "a", params.fuel.a_rf, ...
        "rho_f", params.fuel.rho_f, "web_min", K.C12.lo, "x2", info.x2, ...
        "p_min", params.combustion.p_min, "p_max", params.combustion.p_max, ...
        "time_output", 0:1.0:1.5*K.C6.lo, "fine_ode", params.time.fine_ode);

    S = sizeEngine(info.mdot_ox, ctx, lookupN, cfg);
    if ~S.ok
        error("verifySolver:sizingFailed", "Default design did not size: %s", S.err);
    end

    design = struct("S", S, "ctx", ctx, "lookupN", lookupN, "cfg", cfg);
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
