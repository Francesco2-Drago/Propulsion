function G = compareGeometries(opts)
    % compareGeometries
    % The three-geometry comparison for point (iii) of the assignment: circular
    % port, rod and tube, and star, all at the same port area so that only the
    % TOPOLOGY of the port differs.
    %
    % The point it makes is the perimeter growth law. Under uniform normal
    % regression the perimeter of a port evolves as
    %
    %     dPb/db = 2*pi*(1 - k)
    %
    % with k the number of fuel ISLANDS inside the port. A plain circular port
    % has k = 0 and its perimeter grows at 2*pi (Steiner). A rod and tube has
    % k = 1, one island in the middle: the outer surface grows at 2*pi and the
    % island shrinks at 2*pi, so the total perimeter is very nearly NEUTRAL, and
    % the O/F drift with it. A star has k = 0 but is not convex, so its
    % perimeter grows more slowly than 2*pi while the concave valleys fill in,
    % and then collapses when the tips reach the casing.
    %
    % So: the star buys initial perimeter at the price of losing it during the
    % burn, the rod and tube buys neutrality at the price of a fuel island that
    % has to be held in place, and the circle is what maximizes perimeter GROWTH.
    % INPUT
    %   opts / struct / 1x1   Optional: oxidizer (default "O2(L)"), grid
    %                         (default 350), x1 (default 0.28, the port radius
    %                         ratio the three geometries are matched on),
    %                         N (default 8 star tips), do_plots (default true)
    % OUTPUT
    %   G    / struct / 1x3   One record per geometry: name, k, lookup, Phi0,
    %                         dPb_db, Isp_load, sigma, drift, mdot_ox, Gox0,
    %                         R_c, L

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, "oxidizer"), opts.oxidizer = "O2(L)"; end
    if ~isfield(opts, "grid"), opts.grid = 350; end
    if ~isfield(opts, "x1"), opts.x1 = 0.28; end
    if ~isfield(opts, "N"), opts.N = 8; end
    if ~isfield(opts, "do_plots"), opts.do_plots = true; end

    repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(genpath(repoRoot));

    params = combustion_params();
    params.engine.ext_diameter = 2.0;          % R_c = 1
    params.mdf.grid_divisions = opts.grid;
    params.mdf.perimeter_from_area = false;
    n = params.fuel.n_rf;

    C = optimizationConstraints();
    thermo = get_thermo(opts.oxidizer, repoRoot);
    ctx = engineContext(thermo, params, C);

    % ------------ THE THREE PORTS, MATCHED ON INITIAL AREA ------------
    % Matching the area is what makes the comparison about topology rather than
    % about size: all three start with the same G_ox(0) for a given flow.
    Ap_target = pi * opts.x1^2;                % [-] at R_c = 1

    specs = build_specs(Ap_target, opts.N, params);

    fprintf("\n%s\n", repmat('=', 1, 88));
    fprintf(" GEOMETRY COMPARISON at matched port area Ap(0) = %.5f (x1_eq = %.4f), %s\n", ...
        Ap_target, opts.x1, opts.oxidizer);
    fprintf("%s\n", repmat('-', 1, 88));
    fprintf(" %-12s %3s %8s %8s %9s %9s %8s %8s %8s %8s\n", ...
        "geometry", "k", "Ap(0)", "Pb(0)", "Phi0", "dPb/db", "Isp_load", ...
        "sigma", "drift", "L[m]");
    fprintf("%s\n", repmat('-', 1, 88));

    G = repmat(empty_record(), 1, numel(specs));
    for i = 1:numel(specs)
        G(i) = evaluate_one(specs(i), params, ctx, n);
        g = G(i);
        if isfinite(g.Isp_load)
            fprintf(" %-12s %3d %8.5f %8.5f %9.5f %9.4f %8.2f %8.4f %8.3f %8.3f\n", ...
                g.name, g.k, g.Ap0, g.Pb0, g.Phi0, g.dPb_db, g.Isp_load, ...
                g.sigma, g.drift, g.L);
        else
            fprintf(" %-12s %3d %8.5f %8.5f %9.5f %9.4f %8s %8s %8s %8s   %s\n", ...
                g.name, g.k, g.Ap0, g.Pb0, g.Phi0, g.dPb_db, ...
                "-", "-", "-", "-", g.fail);
        end
    end
    fprintf("%s\n", repmat('-', 1, 88));

    % ------------ THE LAW ------------
    fprintf("\n dPb/db = 2*pi*(1 - k), with k the number of fuel islands in the port,\n");
    fprintf(" fitted over the burnback range where that k actually holds:\n");
    for i = 1:numel(G)
        expected = 2*pi*(1 - G(i).k);
        if isfinite(G(i).k_until)
            range_txt = sprintf("b~ < %.3f, where the island is consumed", G(i).k_until);
        else
            range_txt = sprintf("b~ < %.3f", G(i).b_fit);
        end
        fprintf("   %-12s k = %d  ->  expected %7.4f, measured %7.4f   (%s)\n", ...
            G(i).name, G(i).k, expected, G(i).dPb_db, range_txt);
    end

    fprintf("\n Two things worth reading off this table.\n\n");
    fprintf(" The ROD AND TUBE is neutral only while the rod lasts. The outer\n");
    fprintf(" surface grows at +2*pi and the island shrinks at -2*pi, so they\n");
    fprintf(" cancel; but the rod is thin and burns out early, and from then on\n");
    fprintf(" the port is topologically a circle again and the perimeter resumes\n");
    fprintf(" growing at 2*pi. Neutrality is bought for part of the burn, not all\n");
    fprintf(" of it, and it costs a fuel island that has to be held in place.\n\n");
    fprintf(" The STAR has k = 0 like the circle, so the law predicts 2*pi, but the\n");
    fprintf(" port is NOT CONVEX and the measured slope falls short. Steiner's\n");
    fprintf(" formula holds for convex bodies; a re-entrant valley fills in as it\n");
    fprintf(" regresses, and that perimeter is lost. The star buys initial\n");
    fprintf(" perimeter at the price of losing it during the burn.\n");
    fprintf("%s\n\n", repmat('=', 1, 88));

    if opts.do_plots
        plot_comparison(G, params);
    end
end

%% FUNCTIONS

function e = empty_record()
    % empty_record
    % Prototype of one geometry record.
    % INPUT
    %   None
    % OUTPUT
    %   e / struct / 1x1   Empty record

    e = struct("name", "", "k", 0, "lookup", struct([]), "Ap0", NaN, ...
        "Pb0", NaN, "Phi0", NaN, "dPb_db", NaN, "Isp_load", NaN, ...
        "sigma", NaN, "drift", NaN, "drift_full", NaN, "mdot_ox", NaN, ...
        "Gox0", NaN, "R_c", NaN, "L", NaN, "meshdata", struct([]), ...
        "type", "", "fail", "", "k_until", Inf, "b_fit", NaN);
end

function specs = build_specs(Ap_target, N, params)
    % build_specs
    % The three port geometries, each solved so its initial area matches.
    % INPUT
    %   Ap_target / double / 1x1   Initial port area to match [-] at R_c = 1
    %   N         / double / 1x1   Number of star tips [-]
    %   params    / struct / 1x1   Combustion parameters
    % OUTPUT
    %   specs     / struct / 1x3   name, k, type, meshdata, r_out

    R_c = params.engine.ext_diameter / 2;

    % 1) Circular port: k = 0 fuel islands, convex
    r0 = sqrt(Ap_target/pi);
    s1 = struct("name", "circle", "k", 0, "type", "cylinder", ...
        "meshdata", struct("diameter", 2*r0*R_c), "r_out", r0, ...
        "k_until", Inf);

    % 2) Rod and tube: k = 1, a fuel island in the middle. The annulus between
    %    rod radius ri and tube radius ro must carry the same area.
    %    pi*(ro^2 - ri^2) = Ap_target, with the rod taken at a fixed fraction
    %    of the tube so the geometry is determined.
    %    The island only exists until the rod is consumed, at b = ri: after
    %    that the port is topologically a circle again and k drops to 0. That
    %    is where the neutrality ends, so the slope has to be measured there.
    frac = 0.45;                                    % ri/ro [-]
    ro = sqrt(Ap_target/(pi*(1 - frac^2)));
    ri = frac*ro;
    s2 = struct("name", "rod & tube", "k", 1, "type", "rod_and_tube", ...
        "meshdata", struct("inner_diameter", 2*ri*R_c, ...
                           "outer_diameter", 2*ro*R_c), "r_out", ro, ...
        "k_until", ri);

    % 3) Star: k = 0 but not convex. Area of a star with apothem ri and tip re
    %    is the inner polygon plus N outer triangles; solve for re at a fixed
    %    apothem so the area matches.
    ri_s = 0.80 * r0;                               % [-] a deliberately real star
    A_poly = @(a) N * a^2 * tan(pi/N);
    side = @(a) 2 * a * tan(pi/N);
    A_star = @(re, a) A_poly(a) + N * 0.5 * (re - a) * side(a);
    re_s = fzero(@(re) A_star(re, ri_s) - Ap_target, [ri_s, 4*r0]);
    s3 = struct("name", "star", "k", 0, "type", "star", ...
        "meshdata", struct("inner_diameter", 2*ri_s*R_c, ...
                           "outer_diameter", 2*re_s*R_c, "n_tips", N), ...
        "r_out", re_s, "k_until", Inf);

    specs = [s1, s2, s3];
end

function g = evaluate_one(spec, params, ctx, n)
    % evaluate_one
    % Build one geometry's lookup and score it with the same shapeMerit the
    % optimization uses.
    % INPUT
    %   spec   / struct / 1x1   Geometry specification from build_specs
    %   params / struct / 1x1   Combustion parameters
    %   ctx    / struct / 1x1   Context from engineContext
    %   n      / double / 1x1   Regression exponent [-]
    % OUTPUT
    %   g      / struct / 1x1   Geometry record

    g = empty_record();
    g.name = string(spec.name);
    g.k = spec.k;
    g.type = string(spec.type);
    g.meshdata = spec.meshdata;

    R_c = params.engine.ext_diameter / 2;
    p = params;
    p.geometry.type = spec.type;
    lookup = build_shape_lookup(spec.meshdata, p);
    g.lookup = lookup;

    lookupN.b = lookup.b / R_c;
    lookupN.Ap = lookup.Ap / R_c^2;
    lookupN.perim = lookup.perim / R_c;

    g.Ap0 = lookupN.Ap(1);
    g.Pb0 = lookupN.perim(1);
    g.Phi0 = g.Ap0^n / g.Pb0;

    % Perimeter growth rate, fitted where the declared topology actually holds:
    % away from the wall, where the MDF stops counting surface outside the
    % casing, AND before any fuel island is consumed, after which k drops.
    b = lookupN.b;
    b_fit = min(0.8*b(end), spec.k_until);
    g.k_until = spec.k_until;
    keep = b <= b_fit;
    if nnz(keep) >= 3
        q = polyfit(b(keep), lookupN.perim(keep), 1);
        g.dPb_db = q(1);
    end
    g.b_fit = b_fit;

    merit_cfg = ctx.merit_cfg;
    merit_cfg.r_out = spec.r_out;
    M = shapeMerit(lookupN, ctx.Isp_of, merit_cfg);
    if ~M.ok
        g.fail = "no admissible O/F level (C4/C5)";
        return
    end

    g.Isp_load = M.merit;
    g.sigma = M.sigma;
    g.drift = M.drift;
    g.drift_full = M.drift_full;
    g.mdot_ox = M.mdot_ox;
    g.Gox0 = M.Gox0;
    g.R_c = M.R_c;
    g.L = M.L;
end

function plot_comparison(G, params)
    % plot_comparison
    % Cross-sections, perimeter histories and shape functions of the three
    % ports, side by side.
    % INPUT
    %   G      / struct / 1xM   Geometry records
    %   params / struct / 1x1   Combustion parameters
    % OUTPUT
    %   None (creates a figure)

    R_c = params.engine.ext_diameter / 2;
    n = params.fuel.n_rf;
    nG = numel(G);
    colors = lines(nG);

    figure("Name", "Geometry comparison", "Color", "w");
    tiledlayout(2, nG, "TileSpacing", "compact", "Padding", "compact");

    % Row 1: the cross-sections
    th = linspace(0, 2*pi, 400);
    for i = 1:nG
        nexttile(i)
        fill(R_c*cos(th), R_c*sin(th), [0.82 0.72 0.47], ...
            "EdgeColor", "k", "LineWidth", 1.2);
        hold on
        draw_port(G(i), th);
        axis equal
        axis(R_c*[-1.1 1.1 -1.1 1.1])
        set(gca, "XTick", [], "YTick", [])
        title(sprintf("%s, k = %d", G(i).name, G(i).k));
    end

    % Row 2, left: perimeter against burnback, the growth law
    nexttile(nG + 1, [1 2])
    hold on
    for i = 1:nG
        b = G(i).lookup.b / R_c;
        Pb = G(i).lookup.perim / R_c;
        plot(b, Pb, "LineWidth", 1.6, "Color", colors(i,:), ...
            "DisplayName", sprintf("%s (k=%d, dPb/db=%.2f)", ...
            G(i).name, G(i).k, G(i).dPb_db));
    end
    b_ref = linspace(0, 0.6, 20);
    plot(b_ref, G(1).Pb0 + 2*pi*b_ref, "k--", "LineWidth", 1.1, ...
        "DisplayName", "2\pi slope (k = 0, convex)");
    grid on
    xlabel("$\tilde b$ [-]", "Interpreter", "latex");
    ylabel("$\tilde P_b$ [-]", "Interpreter", "latex");
    title("Perimeter growth: dP_b/db = 2\pi(1 - k)");
    legend("Location", "best", "Interpreter", "none");

    % Row 2, right: the shape function, i.e. the O/F history up to lambda
    nexttile(nG + 3)
    hold on
    for i = 1:nG
        b = G(i).lookup.b / R_c;
        Ap = G(i).lookup.Ap / R_c^2;
        Pb = G(i).lookup.perim / R_c;
        plot(b, Ap.^n ./ max(Pb, realmin), "LineWidth", 1.6, "Color", colors(i,:), ...
            "DisplayName", G(i).name);
    end
    grid on
    xlabel("$\tilde b$ [-]", "Interpreter", "latex");
    ylabel("$\tilde\Phi = \tilde A_p^n/\tilde P_b$ [-]", "Interpreter", "latex");
    title("Shape function (the O/F history)");
    legend("Location", "northwest", "Interpreter", "none");
end

function draw_port(g, th)
    % draw_port
    % Fill one port cross-section, whatever its topology. The casing is already
    % drawn by the caller.
    % INPUT
    %   g   / struct / 1x1   Geometry record
    %   th  / double / 1xM   Angles for the circular outlines [rad]
    % OUTPUT
    %   None (draws into the current axes)

    md = g.meshdata;
    switch g.type
        case "cylinder"
            r = md.diameter/2;
            fill(r*cos(th), r*sin(th), "w", ...
                "EdgeColor", [0.2 0.2 0.2], "LineWidth", 1.1);
        case "rod_and_tube"
            ro = md.outer_diameter/2;
            ri = md.inner_diameter/2;
            fill(ro*cos(th), ro*sin(th), "w", ...
                "EdgeColor", [0.2 0.2 0.2], "LineWidth", 1.1);
            % The rod is the fuel island: same colour as the grain
            fill(ri*cos(th), ri*sin(th), [0.82 0.72 0.47], ...
                "EdgeColor", [0.2 0.2 0.2], "LineWidth", 1.1);
        case "star"
            half = make_mesh0("star", md, 800, "cartesian");
            port = [half; flipud([half(:,1), -half(:,2)])];
            fill(port(:,1), port(:,2), "w", ...
                "EdgeColor", [0.2 0.2 0.2], "LineWidth", 1.1);
        otherwise
            error("draw_port:unknownType", "Unknown geometry type '%s'.", g.type);
    end
end
