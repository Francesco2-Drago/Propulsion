function [J, info] = shapeCostFunction(x1, h, N, params, ctx)
    % shapeCostFunction
    % Phase A objective: MAXIMIZE the specific impulse on the loaded mass of a
    % normalized grain shape at R_c = 1, each shape evaluated at its own best
    % O/F level. Returns J = -merit, or +Inf when a gate fails.
    %
    % The shape is parametrized by (x1, h) instead of (x1, x2): with
    % x2 = x1/cos(pi/N) + h the "valid star" condition C1 is a plain bound on h
    % and the problem is a box again. The old formulation had the optimum
    % sitting exactly on the face of a linear inequality, where direct search
    % converges badly.
    %
    % INPUT
    %   x1     / double / 1x1   Inner apothem ratio ri/R_c [-]
    %   h      / double / 1x1   Tip height above the inner polygon vertex,
    %                           as a fraction of R_c [-]
    %   N      / double / 1x1   Number of star tips (rounded) [-]
    %   params / struct / 1x1   Combustion parameters (from combustion_params())
    %   ctx    / struct / 1x1   Evaluation context built by engineContext:
    %                           K (constraints keyed by id), merit_cfg, Isp_of,
    %                           h_min [-]
    % OUTPUT
    %   J      / double / 1x1   -Isp_load [s], +Inf if any gate fails
    %   info   / struct / 1x1   Diagnostics: x2, merit, Isp_med, OF_med, sigma,
    %                           drift, Gtilde, lambda, Phi0, feasible, fail, err

    info = struct("x2", NaN, "merit", NaN, "Isp_med", NaN, "OF_med", NaN, ...
        "sigma", NaN, "drift", NaN, "Gtilde", NaN, "lambda", NaN, ...
        "Phi0", NaN, "of0", NaN, "of_end", NaN, "b_end", NaN, ...
        "mdot_ox", NaN, "Gox0", NaN, "drift_full", NaN, ...
        "R_c", NaN, "L", NaN, "m_load", NaN, "I_tot", NaN, ...
        "feasible", false, "fail", "", "err", "");
    J = Inf;

    K = ctx.K;
    N = round(N);
    r_vertex = x1 / cos(pi / N);        % [-] inner polygon vertex radius
    x2 = r_vertex + h;                  % [-] tip radius ratio re/R_c
    info.x2 = x2;

    % ------------ C1: VALID STAR ------------
    % Implied by the h >= h_min bound of the search box, checked here so the
    % cost function is safe to call on its own
    if h < K.C1.lo
        info.fail = "C1";
        return
    end

    % ------------ C2: TIP RESOLVABILITY ------------
    % The tip must be at least K.C2.lo grid cells tall, or what gets optimized
    % is discretization noise
    if h < ctx.h_min
        info.fail = "C2";
        return
    end

    % ------------ C3: CONTAINMENT ------------
    if max(x2, r_vertex) > K.C3.hi
        info.fail = "C3";
        return
    end

    try
        % Grain cross-section for the MDF lookup, normalized to R_c = 1
        R_c = params.engine.ext_diameter / 2;     % [m] reference casing radius
        meshdata.inner_diameter = 2 * x1 * R_c;   % [m]
        meshdata.outer_diameter = 2 * x2 * R_c;   % [m]
        meshdata.n_tips = N;

        params.mdf.perimeter_from_area = false;   % contour perimeter
        lookup = build_shape_lookup(meshdata, params);

        % Back to the normalized shape, so the merit is scale-free by construction
        lookupN.b = lookup.b / R_c;
        lookupN.Ap = lookup.Ap / R_c^2;
        lookupN.perim = lookup.perim / R_c;

        % ------------ MERIT, C4, C5 AND C12 ------------
        % shapeMerit enforces the gates from the inside, level by level: C5
        % refuses a shape whose initial O/F cannot reach the CEA domain and
        % ends the burn at the upper crossing, C12 stops the burn on the
        % residual web, and C4 rejects every O/F level whose implied oxidizer
        % flow puts G_ox(0) out of band. A shape is feasible when at least one
        % level survives all of them.
        % The tip is the thinnest point of the web, so it is x2 that sets how
        % much fuel stands between the port and the casing.
        merit_cfg = ctx.merit_cfg;
        merit_cfg.r_out = x2;                     % [-]
        M = shapeMerit(lookupN, ctx.Isp_of, merit_cfg);
        if ~M.ok
            info.fail = "C4/C5";
            return
        end

        info.merit = M.merit;
        info.Isp_med = M.Isp_med;
        info.OF_med = M.OF_med;
        info.sigma = M.sigma;
        info.drift = M.drift;
        info.Gtilde = M.Gtilde;
        info.lambda = M.lambda;
        info.Phi0 = M.Phi0;
        info.of0 = M.of0;
        info.of_end = M.of_end;
        info.b_end = M.b_end;
        info.mdot_ox = M.mdot_ox;
        info.Gox0 = M.Gox0;
        info.drift_full = M.drift_full;
        info.R_c = M.R_c;
        info.L = M.L;
        info.m_load = M.m_load;
        info.I_tot = M.I_tot;

        info.feasible = true;
        info.fail = "OK";
        J = -M.merit;
    catch ME
        % Record the error: an MDF failure and a bad shape must never look alike
        info.err = ME.message;
        info.fail = "error";
        warning("shapeCostFunction:evalFailed", ...
            "Shape evaluation failed at x1 = %.4f, h = %.4f, N = %d: %s", ...
            x1, h, N, ME.message);
        J = Inf;
    end
end
