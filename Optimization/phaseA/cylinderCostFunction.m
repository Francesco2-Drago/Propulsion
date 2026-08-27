function [J, info] = cylinderCostFunction(x1, params, ctx)
    % cylinderCostFunction
    % Phase A objective for a circular port, scored by the SAME shapeMerit as
    % the star family so the two compete on equal terms. One variable, the port
    % radius ratio x1 = r0/R_c.
    %
    % C1 and C2 do not apply: there is no tip to be valid or to be resolvable.
    % C3, C4 and C5 apply exactly as they do to a star.
    %
    % The geometry comes from the SAME MDF lookup the stars use, not from the
    % closed form of cylinderLookup, even though the closed form is exact and a
    % thousand times cheaper. The reason is that the two disagree on I~ by
    % 0.7-1.2 % (grid resolution), I~ enters R_c = (t_b a mdot^n/I~)^(1/(2n+1))
    % and therefore the loaded mass and the objective, and the cylinder-star
    % margin is only 0.2 %. Scoring one family analytically and the other on a
    % grid would make that bias one-sided, and it would be large enough to
    % account for the whole margin on its own. Same discretization for
    % everybody, so whatever error the MDF makes is common mode and cancels in
    % the comparison. cylinderLookup stays as the validation reference.
    % INPUT
    %   x1     / double / 1x1   Port radius ratio r0/R_c [-]
    %   params / struct / 1x1   Combustion parameters (from combustion_params())
    %   ctx    / struct / 1x1   Evaluation context from engineContext
    % OUTPUT
    %   J      / double / 1x1   -Isp_load [s], +Inf if any gate fails
    %   info   / struct / 1x1   Same diagnostics as shapeCostFunction

    info = struct("x2", NaN, "merit", NaN, "Isp_med", NaN, "OF_med", NaN, ...
        "sigma", NaN, "drift", NaN, "Gtilde", NaN, "lambda", NaN, ...
        "Phi0", NaN, "of0", NaN, "of_end", NaN, "b_end", NaN, ...
        "mdot_ox", NaN, "Gox0", NaN, "drift_full", NaN, ...
        "R_c", NaN, "L", NaN, "m_load", NaN, "I_tot", NaN, ...
        "feasible", false, "fail", "", "err", "");
    J = Inf;

    % A circular port has no tip: its outer radius is the port radius itself
    info.x2 = x1;

    % ------------ C3: CONTAINMENT ------------
    if x1 <= 0 || x1 >= ctx.K.C3.hi
        info.fail = "C3";
        return
    end

    try
        R_c = params.engine.ext_diameter / 2;     % [m] reference casing radius

        % Same MDF path as the star family, see the note above
        params_c = params;
        params_c.geometry.type = "cylinder";
        params_c.mdf.perimeter_from_area = false;   % contour perimeter
        meshdata.diameter = 2 * x1 * R_c;           % [m]
        lookup = build_shape_lookup(meshdata, params_c);

        % Back to the normalized shape, so the merit is scale-free
        lookupN.b = lookup.b / R_c;
        lookupN.Ap = lookup.Ap / R_c^2;
        lookupN.perim = lookup.perim / R_c;

        % ------------ MERIT, C4, C5 AND C12 ------------
        merit_cfg = ctx.merit_cfg;
        merit_cfg.r_out = x1;                     % [-] the port is the whole grain
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
        info.drift_full = M.drift_full;
        info.Gtilde = M.Gtilde;
        info.lambda = M.lambda;
        info.Phi0 = M.Phi0;
        info.of0 = M.of0;
        info.of_end = M.of_end;
        info.b_end = M.b_end;
        info.mdot_ox = M.mdot_ox;
        info.Gox0 = M.Gox0;
        info.R_c = M.R_c;
        info.L = M.L;
        info.m_load = M.m_load;
        info.I_tot = M.I_tot;

        info.feasible = true;
        info.fail = "OK";
        J = -M.merit;
    catch ME
        info.err = ME.message;
        info.fail = "error";
        warning("cylinderCostFunction:evalFailed", ...
            "Cylinder evaluation failed at x1 = %.4f: %s", x1, ME.message);
        J = Inf;
    end
end
