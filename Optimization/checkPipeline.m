function ok = checkPipeline(opts)
    % checkPipeline
    % Fast smoke test of the whole optimization chain. Runs in a couple of
    % minutes instead of the hour and a half a full main_optimization takes,
    % and is meant to answer one question: is anything obviously broken?
    %
    % It exercises every module on a reduced problem (few tip counts, coarse
    % grid, one oxidizer for the end-to-end part) and checks the analytic
    % identities of the brief, which are what actually catch physics errors.
    % A green run does NOT prove the optimum is right; it proves the machinery
    % is wired correctly and the physics closes.
    % INPUT
    %   opts / struct / 1x1   Optional: oxidizer (default "O2(L)"), grid
    %                         (default 350), full (default false, run the
    %                         reduced phase A + phase B end to end)
    % OUTPUT
    %   ok   / logical / 1x1  True when every check passed

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, "oxidizer"), opts.oxidizer = "O2(L)"; end
    if ~isfield(opts, "grid"), opts.grid = 350; end
    if ~isfield(opts, "full"), opts.full = false; end

    % This file sits at Optimization/, so the root is two levels up
    repoRoot = fileparts(fileparts(mfilename("fullpath")));
    addpath(genpath(repoRoot));

    T = tic;
    R = results_new();

    fprintf("\n%s\n", repmat('=', 1, 72));
    fprintf(" PIPELINE CHECK - oxidizer %s, grid %d\n", opts.oxidizer, opts.grid);
    fprintf("%s\n", repmat('=', 1, 72));

    params = combustion_params();
    params.geometry.type = "star";
    params.engine.ext_diameter = 2.0;          % R_c = 1
    params.mdf.grid_divisions = opts.grid;
    params.mdf.perimeter_from_area = false;
    n = params.fuel.n_rf;
    a = params.fuel.a_rf;

    % ------------ 1. CONSTRAINT TABLE ------------
    fprintf("\n--- 1. Constraint table ---\n");
    [C, K, box] = optimizationConstraints();
    R = check(R, "C1..C14 present and in order", ...
        numel(C) == 14 && all([C.id] == "C" + string(1:14)), ...
        sprintf("%d rows", numel(C)));
    R = check(R, "C6 = 300 s, C7 = 50 kN, C8 = 20 bar", ...
        K.C6.lo == 300 && K.C7.lo == 50e3 && K.C8.lo == 20e5, "");
    R = check(R, "C4 is C10 with a 5% margin", ...
        abs(K.C4.lo/K.C10.lo - 1.05) < 1e-9 && abs(K.C10.hi/K.C4.hi - 1.05) < 1e-9, ...
        sprintf("C4 = [%.0f, %.0f]", K.C4.lo, K.C4.hi));
    R = check(R, "C14 withdrawn, oxidizer free", ...
        K.C14.kind == "removed", "");
    R = check(R, "no mdot_ox search box left", ...
        ~isfield(box, "mdot_ox"), "");
    R = check(R, "every constraint declares who consumes it", ...
        all(arrayfun(@(c) strlength(c.interaction) > 0, C)), ...
        sprintf("%d of %d", sum(arrayfun(@(c) strlength(c.interaction) > 0, C)), numel(C)));

    % C2's floor is grid-dependent, and phase A uses two grids. The star sits
    % exactly on whichever floor is in force while the cylinder has no such
    % bound, so a mismatch between the search grid and the ranking grid biases
    % the comparison towards the cylinder. Assert the floors are what they
    % should be at both resolutions.
    R = check(R, "C2 floor halves from the search grid to the ranking grid", ...
        abs((K.C2.lo/350)/(K.C2.lo/700) - 2) < 1e-12, ...
        sprintf("%.5f at 350, %.5f at 700", K.C2.lo/350, K.C2.lo/700));

    % C4 at 10 kg/s must be the Gtilde window the brief quotes
    Cs = 10^(1/(2*n+1)) / (K.C6.lo*a)^(2/(2*n+1));
    R = check(R, "at mdot_ox = 10, C10 means Gtilde in [3.38, 5.91]", ...
        abs(K.C10.lo/Cs - 3.38) < 0.01 && abs(K.C10.hi/Cs - 5.91) < 0.01, ...
        sprintf("[%.2f, %.2f]", K.C10.lo/Cs, K.C10.hi/Cs));

    % ------------ 2. THERMOCHEMISTRY AND NOZZLE ------------
    fprintf("\n--- 2. Thermochemistry and nozzle ---\n");
    thermo = get_thermo(opts.oxidizer, repoRoot);
    ctx = engineContext(thermo, params, C);
    P = performancePoint(2.0, 20, ctx.thermo, ctx.engine);
    R = check(R, "performancePoint sane at O/F = 2, 20 bar", ...
        P.Tc > 2000 && P.Tc < 4500 && P.Me > 3 && P.Me < 8 && ...
        P.cstar > 1200 && P.cstar < 2200 && P.Isp > 250 && P.Isp < 400, ...
        sprintf("Tc %.0f K, Me %.2f, cstar %.0f m/s, Isp %.1f s", ...
        P.Tc, P.Me, P.cstar, P.Isp));
    R = check(R, "Isp = c_eff/g0 consistent with ceff_of", ...
        abs(ctx.ceff_of(2.0)/(9.80665*ctx.Isp_of(2.0)) - 1) < 1e-6, "");
    R = check(R, "Isp peak inside the CEA domain", ...
        ctx.of_best > K.C5.lo && ctx.of_best < K.C5.hi, ...
        sprintf("peak %.2f s at O/F = %.2f", ctx.Isp_of(ctx.of_best), ctx.of_best));

    % ------------ 3. GEOMETRY: ANALYTIC vs MDF ------------
    fprintf("\n--- 3. Geometry, analytic cylinder vs MDF ---\n");
    x1 = 0.28;
    la = cylinderLookup(x1, params);
    pc = params; pc.geometry.type = "cylinder";
    lm = build_shape_lookup(struct("diameter", 2*x1), pc);
    R = check(R, "Ap(0) agrees within 1%", ...
        rel(lm.Ap(1), la.Ap(1)) < 0.01, sprintf("%.2f%%", 100*rel(lm.Ap(1), la.Ap(1))));
    R = check(R, "Pb(0) agrees within 1%", ...
        rel(lm.perim(1), la.perim(1)) < 0.01, ...
        sprintf("%.2f%%", 100*rel(lm.perim(1), la.perim(1))));
    Ia = trapz(la.b, la.Ap.^n);
    Im = trapz(lm.b, lm.Ap.^n);
    R = check(R, "burn integral I~ agrees within 2%", rel(Im, Ia) < 0.02, ...
        sprintf("%.2f%%", 100*rel(Im, Ia)));

    % Steiner, on the useful part: the MDF zeroes the perimeter at the wall
    fr = lm.b / lm.b(end);
    m = fr <= 0.98;
    p = polyfit(lm.b(m), lm.perim(m), 1);
    R = check(R, "Steiner dPb/db = 2*pi within 2% (useful part)", ...
        rel(p(1), 2*pi) < 0.02, sprintf("%.5f vs %.5f", p(1), 2*pi));

    % ------------ 4. SCALE-FREE IDENTITIES ------------
    fprintf("\n--- 4. Scale-free identities ---\n");
    % Gtilde of a circular port, closed form against the brief's number
    r0 = 0.2638;
    lc = cylinderLookup(r0, params);
    I_c = trapz(lc.b, lc.Ap.^n);
    Gt = I_c^(2/(2*n+1)) / lc.Ap(1);
    R = check(R, "Gtilde of the r0/Rc = 0.2638 port is 4.242", ...
        abs(Gt - 4.242) < 0.01, sprintf("%.4f", Gt));

    % Cylindrical drift, (Rc/r0)^(2n-1)
    Phi_c = lc.Ap.^n ./ lc.perim;
    R = check(R, "cylindrical drift = (Rc/r0)^(2n-1) within 1%", ...
        rel(Phi_c(end)/Phi_c(1), (1/r0)^(2*n-1)) < 0.01, ...
        sprintf("%.4f vs %.4f", Phi_c(end)/Phi_c(1), (1/r0)^(2*n-1)));

    % ------------ 5. PHASE A EVALUATORS ------------
    fprintf("\n--- 5. Phase A evaluators ---\n");
    [Jc, ic] = cylinderCostFunction(0.28, params, ctx);
    R = check(R, "cylinderCostFunction returns a feasible design", ...
        isfinite(Jc) && ic.feasible, sprintf("Isp_load %.2f s, %s", -Jc, ic.fail));

    h = K.C2.lo / opts.grid;
    [Js, is] = shapeCostFunction(0.27, h, 18, params, ctx);
    R = check(R, "shapeCostFunction returns a feasible design", ...
        isfinite(Js) && is.feasible, sprintf("Isp_load %.2f s, %s", -Js, is.fail));

    if isfinite(Jc) && isfinite(Js)
        % The decomposition must reproduce the direct objective
        for nm = ["cylinder", "star"]
            if nm == "cylinder", i = ic; else, i = is; end
            rhs = i.Isp_med*(i.OF_med + 1)/(i.OF_med + 1/(1 - i.sigma));
            R = check(R, "Isp_load decomposition, " + nm + ", within 0.1%", ...
                rel(i.merit, rhs) < 1e-3, sprintf("%.4f vs %.4f", i.merit, rhs));
        end

        % Gox0 must be the 2.4 factorization evaluated on the computed flow
        Cs2 = ic.mdot_ox^(1/(2*n+1)) / (K.C6.lo*a)^(2/(2*n+1));
        R = check(R, "Gox0 = C(mdot_ox)*Gtilde within 1%", ...
            rel(ic.Gox0, Cs2*ic.Gtilde) < 0.01, ...
            sprintf("%.1f vs %.1f", ic.Gox0, Cs2*ic.Gtilde));

        % C6 must hold on the phase A prediction. I~ is recovered from Gtilde
        % and Ap~(0), so this closes the loop through a different route than
        % the one that produced R_c.
        Ap0_n = pi * ic.x2^2;                     % circular port at R_c = 1
        I_rec = (ic.Gtilde * Ap0_n)^((2*n + 1)/2);
        R_c_exp = (K.C6.lo * a * ic.mdot_ox^n / I_rec)^(1/(2*n + 1));
        R = check(R, "R_c satisfies the burn-time equation C6 within 1%", ...
            rel(ic.R_c, R_c_exp) < 0.01, ...
            sprintf("%.4f m vs %.4f m", ic.R_c, R_c_exp));

        % Both gates must actually be inside their band
        R = check(R, "phase A G_ox(0) inside C4", ...
            ic.Gox0 >= K.C4.lo && ic.Gox0 <= K.C4.hi, sprintf("%.0f", ic.Gox0));

        % Perimeter collapse at the casing. The MDF drives Pb to zero at the
        % wall for EVERY shape, cylinder included, so drift_full diverges for
        % both and its value alone says nothing. What separates the families is
        % how WIDE the collapse zone is: for a circular port every point of the
        % surface reaches the casing at the same instant, so it is a couple of
        % lookup levels; for a star the tips arrive well before the flats.
        % What actually separates the two families is the SLIVER at the common
        % stop criterion C12: the star's tips reach the casing before its
        % flats, so when the thinnest point of the web is down to 3 mm there is
        % still fuel left elsewhere. This is the term that carries the margin.
        R = check(R, "star leaves more sliver than the cylinder at the C12 stop", ...
            is.sigma > ic.sigma, ...
            sprintf("star %.4f, cylinder %.4f, gap %.4f", ...
            is.sigma, ic.sigma, is.sigma - ic.sigma));

        % Reported, not gated: the perimeter shrinks near the wall for BOTH
        % families, because the MDF stops counting surface outside the casing.
        % The zones are comparable, so this is not what decides the comparison.
        w_cyl = collapse_width(0.28, params, "cylinder", 0, 0);
        w_star = collapse_width(0.27, params, "star", h, 18);
        fprintf("  [info] %-52s star %.2f%%, cylinder %.2f%% of the web\n", ...
            "perimeter shrinking zone", 100*w_star, 100*w_cyl);

        % Sanity of the predicted engine
        R = check(R, "phase A predicts a sane engine", ...
            2*ic.R_c > 0.4 && 2*ic.R_c < 1.2 && ic.L > 2 && ic.L < 12, ...
            sprintf("2R_c %.3f m, L %.3f m, m_load %.0f kg", ...
            2*ic.R_c, ic.L, ic.m_load));
    end

    % ------------ 6. PHASE B SIZING ------------
    fprintf("\n--- 6. Phase B sizing (one flow, exact ODE) ---\n");
    if isfinite(Jc)
        pn = params; pn.geometry.type = "cylinder";
        lookupN = build_shape_lookup(struct("diameter", 2*ic.x2), pn);
        cfg = struct("n", n, "a", a, "rho_f", params.fuel.rho_f, ...
            "web_min", K.C12.lo, "x2", ic.x2, ...
            "p_min", params.combustion.p_min, "p_max", params.combustion.p_max, ...
            "time_output", 0:1.0:1.5*K.C6.lo, "fine_ode", params.time.fine_ode);
        S = sizeEngine(ic.mdot_ox, ctx, lookupN, cfg);
        R = check(R, "sizeEngine converges", S.ok, string(S.err));
        if S.ok
            R = check(R, "C6 met: burn time = 300 s within 1%", ...
                rel(S.burn_time, K.C6.lo) < 0.01, sprintf("%.1f s", S.burn_time));
            R = check(R, "C7 met: thrust = 50 kN within 1%", ...
                rel(S.mean_thrust, K.C7.lo) < 0.01, ...
                sprintf("%.2f kN mean, %.2f kN initial", ...
                S.mean_thrust*1e-3, S.thrust0*1e-3));
            R = check(R, "C8 met: peak p_c = 20 bar within 1%", ...
                rel(S.p_peak, K.C8.lo) < 0.01, sprintf("%.2f bar", S.p_peak*1e-5));
            R = check(R, "C12 met: residual web >= 3 mm", ...
                S.web_residual >= K.C12.lo - 1e-9, ...
                sprintf("%.1f mm", S.web_residual*1e3));

            % Mass conservation identity of the brief 2.8
            lhs = 1/S.Gox_end;
            rhs = 1/S.Gox0 + S.burn_time/(params.fuel.rho_f*S.L*S.OF_med);
            R = check(R, "1/Gox_end = 1/Gox0 + t_b/(rho_f L OF_med) within 2%", ...
                rel(lhs, rhs) < 0.02, sprintf("%.5f vs %.5f", lhs, rhs));

            % Phase A predicted the size before the ODE ever ran
            R = check(R, "phase A R_c predicted phase B within 3%", ...
                rel(S.R_c, ic.R_c) < 0.03, ...
                sprintf("A %.4f m, B %.4f m", ic.R_c, S.R_c));
            R = check(R, "phase A L predicted phase B within 8%", ...
                rel(S.L, ic.L) < 0.08, sprintf("A %.3f m, B %.3f m", ic.L, S.L));
            R = check(R, "phase A Isp_load predicted phase B within 3%", ...
                rel(S.Isp_load, -Jc) < 0.03, ...
                sprintf("A %.2f s, B %.2f s", -Jc, S.Isp_load));
        end
    end

    % ------------ 7. END TO END, REDUCED ------------
    if opts.full
        fprintf("\n--- 7. End to end on a reduced problem ---\n");
        optsA = struct("do_plots", false, "quiet", true, "write_csv", false, ...
            "n_starts", 1, "n_report", 3, "grid_search", opts.grid, ...
            "grid_fine", opts.grid);
        shapes = phaseA(thermo, params, C, optsA);
        R = check(R, "phaseA returns a ranking", ~isempty(shapes), ...
            sprintf("%d shapes, best %s at %.2f s", numel(shapes), ...
            best_geom(shapes), best_isp(shapes)));
        if ~isempty(shapes)
            design = phaseB(shapes, thermo, params, C, ...
                struct("do_plots", false, "quiet", true, "n_shapes", 2, ...
                "n_sweep", 3, "validate", false));
            R = check(R, "phaseB returns a design", design.ok, string(design.failTag));
            if design.ok
                R = check(R, "end-to-end design is sane", ...
                    2*design.S.R_c > 0.4 && 2*design.S.R_c < 1.2 && ...
                    design.S.L > 2 && design.S.L < 12 && design.S.Isp_load > 300, ...
                    sprintf("%s, 2R_c %.3f m, L %.3f m, Isp_load %.2f s", ...
                    design.shape.geometry, 2*design.S.R_c, design.S.L, ...
                    design.S.Isp_load));
            end
        end
    else
        fprintf("\n--- 7. End to end: SKIPPED (pass opts.full = true) ---\n");
    end

    % ------------ SUMMARY ------------
    ok = R.failed == 0;
    fprintf("\n%s\n", repmat('=', 1, 72));
    if ok
        fprintf(" ALL %d CHECKS PASSED in %.0f s\n", R.total, toc(T));
    else
        fprintf(" %d of %d CHECKS FAILED in %.0f s\n", R.failed, R.total, toc(T));
        for i = 1:numel(R.failures)
            fprintf("   FAIL: %s\n", R.failures(i));
        end
    end
    fprintf("%s\n\n", repmat('=', 1, 72));
end

%% FUNCTIONS

function R = results_new()
    % results_new
    % Empty tally of the checks.
    % INPUT
    %   None
    % OUTPUT
    %   R / struct / 1x1   total, failed, failures

    R = struct("total", 0, "failed", 0, "failures", strings(1, 0));
end

function R = check(R, name, condition, detail)
    % check
    % Record one check and print it.
    % INPUT
    %   R         / struct  / 1x1   Tally
    %   name      / string  / 1x1   What is being checked
    %   condition / logical / 1x1   Outcome
    %   detail    / string  / 1x1   Measured value, printed either way
    % OUTPUT
    %   R         / struct  / 1x1   Updated tally

    R.total = R.total + 1;
    if condition
        tag = "ok  ";
    else
        tag = "FAIL";
        R.failed = R.failed + 1;
        R.failures(end+1) = name + "  (" + string(detail) + ")";
    end
    if strlength(string(detail)) > 0
        fprintf("  [%s] %-52s %s\n", tag, name, detail);
    else
        fprintf("  [%s] %s\n", tag, name);
    end
end

function d = rel(a, b)
    % rel
    % Relative difference, safe at zero.
    % INPUT
    %   a, b / double / 1x1   Values
    % OUTPUT
    %   d    / double / 1x1   |a/b - 1|, Inf when b is zero and a is not

    if b == 0
        d = double(a ~= 0) * Inf;
    else
        d = abs(a/b - 1);
    end
end

function w = collapse_width(x1, params, geometry, h, N)
    % collapse_width
    % Fraction of the web over which the burning perimeter is SHRINKING. Under
    % uniform normal regression the perimeter of a port can only grow, so
    % dPb/db < 0 means the surface is being eaten by the casing: this isolates
    % the wall collapse from the slower-than-Steiner growth that any non-convex
    % port shows anyway while its concave corners fill in.
    %
    % The MDF drives Pb to zero at the wall for every shape, cylinder included,
    % so what separates the families is how early the shrinking starts.
    % INPUT
    %   x1       / double / 1x1   Port ratio [-]
    %   params   / struct / 1x1   Combustion parameters
    %   geometry / string / 1x1   "cylinder" or "star"
    %   h, N     / double / 1x1   Tip height and count, stars only
    % OUTPUT
    %   w        / double / 1x1   Collapse width, as a fraction of the web [-]

    p = params;
    if geometry == "cylinder"
        p.geometry.type = "cylinder";
        lk = build_shape_lookup(struct("diameter", 2*x1), p);
    else
        p.geometry.type = "star";
        x2 = x1/cos(pi/N) + h;
        lk = build_shape_lookup(struct("inner_diameter", 2*x1, ...
            "outer_diameter", 2*x2, "n_tips", N), p);
    end

    b = lk.b(:);
    Pb = lk.perim(:);
    dPb = gradient(Pb, b);
    k = find(dPb < 0, 1, "first");
    if isempty(k)
        w = 0;
    else
        w = (b(end) - b(k)) / b(end);
    end
end

function g = best_geom(shapes)
    % best_geom
    % Geometry of the best shape in a ranking.
    % INPUT
    %   shapes / struct / 1xK   Ranking
    % OUTPUT
    %   g      / string / 1x1   Geometry name

    if isempty(shapes)
        g = "-";
    else
        g = string(shapes(1).geometry);
    end
end

function v = best_isp(shapes)
    % best_isp
    % Objective of the best shape in a ranking.
    % INPUT
    %   shapes / struct / 1xK   Ranking
    % OUTPUT
    %   v      / double / 1x1   Isp on the loaded mass [s]

    if isempty(shapes)
        v = NaN;
    else
        v = shapes(1).Isp_load;
    end
end
