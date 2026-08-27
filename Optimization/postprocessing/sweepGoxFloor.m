function W = sweepGoxFloor(opts)
    % sweepGoxFloor
    % The curve of the brief 6.2: the objective against the initial oxidizer
    % flux, swept from below the C10 floor up to its ceiling.
    %
    % Why it exists. C10 is the constraint that closes the problem, and its
    % floor ALWAYS binds: without it the optimum runs away towards a thin fuel
    % annulus in an enormous casing. A binding constraint whose value was picked
    % by hand has to be justified, and the honest way to justify it is to show
    % what the design does on either side of it, not to bury a weight in the
    % cost function. That is the whole argument for having no weighted penalties
    % anywhere in this optimization.
    %
    % The sweep deliberately goes BELOW the C10 floor, so it has to relax C4 to
    % get there. That is legitimate here and only here: this is a diagnostic
    % study, not the optimization. C10 itself is untouched.
    % INPUT
    %   opts / struct / 1x1   Optional: oxidizer (default "O2(L)"), grid
    %                         (default 350), Gox_range (default [200 700]),
    %                         n_points (default 11), do_plots (default true)
    % OUTPUT
    %   W    / struct / 1x1   Gox0, Isp_load, R_c, L, sigma, drift, OF_med,
    %                         mdot_ox, x1, all 1xM, plus the C10 band

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, "oxidizer"), opts.oxidizer = "O2(L)"; end
    if ~isfield(opts, "grid"), opts.grid = 350; end
    if ~isfield(opts, "Gox_range"), opts.Gox_range = [200, 700]; end
    if ~isfield(opts, "n_points"), opts.n_points = 11; end
    if ~isfield(opts, "do_plots"), opts.do_plots = true; end

    repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(genpath(repoRoot));

    params = combustion_params();
    params.geometry.type = "cylinder";
    params.engine.ext_diameter = 2.0;          % R_c = 1
    params.mdf.grid_divisions = opts.grid;
    params.mdf.perimeter_from_area = false;

    C = optimizationConstraints();
    K = constraintsById(C);
    thermo = get_thermo(opts.oxidizer, repoRoot);
    ctx = engineContext(thermo, params, C);

    % Relax C4 so the sweep can reach below the C10 floor. Only the diagnostic
    % context is relaxed; the constraint table is not touched.
    ctx_free = ctx;
    ctx_free.merit_cfg.Gox_lo = 0;
    ctx_free.merit_cfg.Gox_hi = Inf;

    % G_ox(0) is monotone in the port radius, so sweeping the port IS sweeping
    % the flux: no inversion needed. Doing it the other way round, bisecting on
    % x1 for each flux target, would cost forty MDF solves per point.
    x_scan = linspace(0.10, 0.45, max(2*opts.n_points, 20));
    m = numel(x_scan);
    [Gox0, Isp_load, R_c, L, sigma, drift, OF_med, mdot_ox, x1] = deal(nan(1, m));

    fprintf("\n%s\n", repmat('=', 1, 92));
    fprintf(" Isp ON LOADED MASS vs INITIAL OXIDIZER FLUX - %s, circular port\n", ...
        opts.oxidizer);
    fprintf(" C10 band [%g, %g] kg/(m2 s); the sweep crosses the floor on purpose\n", ...
        K.C10.lo, K.C10.hi);
    fprintf("%s\n", repmat('-', 1, 92));
    fprintf(" %8s %9s %9s %8s %8s %8s %8s %8s %8s\n", ...
        "x1", "Gox0", "Isp_load", "2R_c[m]", "L[m]", "L/D", "sigma", ...
        "drift", "O/F med");
    fprintf("%s\n", repmat('-', 1, 92));

    for i = 1:m
        [~, info] = cylinderCostFunction(x_scan(i), params, ctx_free);
        if ~info.feasible
            continue
        end
        x1(i) = x_scan(i);
        Gox0(i) = info.Gox0;
        Isp_load(i) = info.merit;
        R_c(i) = info.R_c;
        L(i) = info.L;
        sigma(i) = info.sigma;
        drift(i) = info.drift;
        OF_med(i) = info.OF_med;
        mdot_ox(i) = info.mdot_ox;
    end

    % Keep only what falls in the requested flux window, and order by flux
    inband = isfinite(Gox0) & Gox0 >= opts.Gox_range(1) & Gox0 <= opts.Gox_range(2);
    [Gox0, idx] = sort(Gox0(inband));
    take = @(v) v(inband);
    Isp_load = subs_sorted(take(Isp_load), idx);
    R_c = subs_sorted(take(R_c), idx);
    L = subs_sorted(take(L), idx);
    sigma = subs_sorted(take(sigma), idx);
    drift = subs_sorted(take(drift), idx);
    OF_med = subs_sorted(take(OF_med), idx);
    mdot_ox = subs_sorted(take(mdot_ox), idx);
    x1 = subs_sorted(take(x1), idx);

    for i = 1:numel(Gox0)
        fprintf(" %8.4f %9.0f %9.2f %8.3f %8.3f %8.2f %8.4f %8.3f %8.3f\n", ...
            x1(i), Gox0(i), Isp_load(i), 2*R_c(i), L(i), L(i)/(2*R_c(i)), ...
            sigma(i), drift(i), OF_med(i));
    end
    fprintf("%s\n", repmat('-', 1, 92));

    W = struct("Gox0", Gox0, "Isp_load", Isp_load, "R_c", R_c, "L", L, ...
        "sigma", sigma, "drift", drift, "OF_med", OF_med, ...
        "mdot_ox", mdot_ox, "x1", x1, ...
        "C10", [K.C10.lo, K.C10.hi], "oxidizer", string(opts.oxidizer));

    % What the curve says, in numbers
    good = isfinite(Isp_load);
    if nnz(good) >= 2
        [~, i_top] = max(Isp_load);
        fprintf("\n Best Isp_load %.2f s at G_ox(0) = %.0f", Isp_load(i_top), Gox0(i_top));
        if Gox0(i_top) < K.C10.lo
            fprintf(", i.e. BELOW the C10 floor: the floor costs %.2f s\n", ...
                Isp_load(i_top) - interp_at(Gox0, Isp_load, K.C10.lo));
        else
            fprintf("\n");
        end
        fprintf(" Over the swept range the casing moves %.3f -> %.3f m and L %.3f -> %.3f m,\n", ...
            2*min(R_c(good)), 2*max(R_c(good)), min(L(good)), max(L(good)));
        fprintf(" while the drift moves %.2f -> %.2f: above ~200 the geometry is nearly\n", ...
            min(drift(good)), max(drift(good)));
        fprintf(" insensitive and what really changes is the mixture-ratio excursion.\n");
    end
    fprintf("%s\n\n", repmat('=', 1, 92));

    if opts.do_plots
        plot_sweep(W, K);
    end
end

%% FUNCTIONS

function v = subs_sorted(v, idx)
    % subs_sorted
    % Reorder a filtered sweep column the same way the flux column was sorted.
    % INPUT
    %   v   / double / 1xM   Column already restricted to the in-band points
    %   idx / double / 1xM   Sort permutation of the flux column
    % OUTPUT
    %   v   / double / 1xM   Reordered column

    v = v(idx);
end

function y = interp_at(x, y_vals, x0)
    % interp_at
    % Linear interpolation of a sweep at one abscissa, ignoring gaps.
    % INPUT
    %   x      / double / 1xM   Abscissae
    %   y_vals / double / 1xM   Ordinates
    %   x0     / double / 1x1   Where to interpolate
    % OUTPUT
    %   y      / double / 1x1   Interpolated value, NaN if not possible

    good = isfinite(x) & isfinite(y_vals);
    if nnz(good) < 2
        y = NaN;
        return
    end
    y = interp1(x(good), y_vals(good), x0, "linear", "extrap");
end

function plot_sweep(W, K)
    % plot_sweep
    % The figure of the brief 6.2: the objective and the geometry against the
    % initial flux, with the C10 band shaded.
    % INPUT
    %   W / struct / 1x1   Sweep results
    %   K / struct / 1x1   Constraints keyed by id
    % OUTPUT
    %   None (creates a figure)

    good = isfinite(W.Isp_load);
    if nnz(good) < 2
        return
    end

    figure("Name", "Isp on loaded mass vs initial G_ox", "Color", "w");
    tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");

    nexttile
    shade_band(K.C10.lo, K.C10.hi);
    hold on
    plot(W.Gox0(good), W.Isp_load(good), "o-", "LineWidth", 1.6, ...
        "Color", [0 0.45 0.74], "DisplayName", "I_{sp} on loaded mass");
    grid on
    xlabel("G_{ox}(0) [kg/(m^2 s)]");
    ylabel("I_{sp} on loaded mass [s]");
    title("The objective");
    legend("Location", "best");

    nexttile
    shade_band(K.C10.lo, K.C10.hi);
    hold on
    yyaxis left
    plot(W.Gox0(good), 2*W.R_c(good), "o-", "LineWidth", 1.5);
    ylabel("casing diameter 2R_c [m]");
    yyaxis right
    plot(W.Gox0(good), W.L(good), "s-", "LineWidth", 1.5);
    ylabel("grain length L [m]");
    grid on
    xlabel("G_{ox}(0) [kg/(m^2 s)]");
    title("The geometry barely moves");

    nexttile
    shade_band(K.C10.lo, K.C10.hi);
    hold on
    plot(W.Gox0(good), W.drift(good), "o-", "LineWidth", 1.6, ...
        "Color", [0.85 0.33 0.10]);
    grid on
    xlabel("G_{ox}(0) [kg/(m^2 s)]");
    ylabel("O/F drift over the useful burn [-]");
    title("What actually changes");
end

function shade_band(lo, hi)
    % shade_band
    % Shade the admissible C10 band behind a plot.
    % INPUT
    %   lo, hi / double / 1x1   Band edges [kg/(m2 s)]
    % OUTPUT
    %   None (draws into the current axes)

    yl = [-1e6, 1e6];
    patch([lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], [0.90 0.94 0.90], ...
        "EdgeColor", "none", "HandleVisibility", "off");
    hold on
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
