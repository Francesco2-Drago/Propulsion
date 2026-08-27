function shapes = phaseA(thermo, params, C, opts)
% phaseA
% PHASE A - normalized grain shape optimization at R_c = 1.
%
% Maximize the specific impulse on the LOADED mass over the normalized star
% family, each shape evaluated at its own best O/F level. Nothing here knows
% the size of the engine: everything computable at R_c = 1 lives in this phase,
% including the initial-flux gate C4, which is C10 evaluated on the oxidizer
% flow the thrust requirement implies.
%
% The integer tip count N is enumerated in parallel; for each N the continuous
% problem (x1, h) = (ri/R_c, tip height) is solved with patternsearch from
% several starting points. The output is a RANKING, not a single winner: phase
% B tries the shapes in order until one passes its own gates.
%
% The thermochemistry is an INPUT: the oxidizer parametrizes the whole A -> B
% chain, because it moves the peak of Isp(O/F), hence the best O/F level, hence
% mdot_ox, hence G_ox(0), hence the admissible window on the port fraction, and
% so potentially the optimal shape itself.
% INPUT
%   thermo / struct / 1x1   Oxidizer thermochemistry, from get_thermo
%   params / struct / 1x1   Combustion parameters (from combustion_params())
%   C      / struct / 1xM   Constraint table from optimizationConstraints
%   opts   / struct / 1x1   Optional: do_plots (default false), quiet
%                           (default false), write_csv (default true),
%                           n_starts, n_report, grid_search, grid_fine
% OUTPUT
%   shapes / struct / 1xK   Ranked shapes, best first. Fields: N, x1, h, x2,
%                           Isp_load, Isp_med, OF_med, sigma, drift, Gtilde,
%                           lambda, Phi0, of0, of_end, b_end, mdot_ox, Gox0

if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = default_opt(opts, "do_plots", false);
opts = default_opt(opts, "quiet", false);
opts = default_opt(opts, "write_csv", true);
opts = default_opt(opts, "n_starts", 3);          % [-] starting points per tip count
opts = default_opt(opts, "n_report", 10);         % [-] rows of the final ranking
opts = default_opt(opts, "grid_search", 350);     % [-] coarse grid for the search
opts = default_opt(opts, "grid_fine", 700);       % [-] grid for the re-evaluation

repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));

params.geometry.type = "star";
params.engine.ext_diameter = 2.0;    % [m] normalized casing, R_c = 1
params.mdf.grid_divisions = opts.grid_search;
params.mdf.perimeter_from_area = false;

n_starts = opts.n_starts;
n_report = opts.n_report;
grid_fine = opts.grid_fine;
resultFile = fullfile(repoRoot, "shapeOptimizationResults.csv");

% Shared context: constraints, thermochemistry, Isp and c_eff curves
ctx = engineContext(thermo, params, C);

% Search box 
lb = [ctx.box.x1(1); ctx.h_min];
ub = [ctx.box.x1(2); 1 - ctx.box.x1(1)];
N_list = ctx.box.N(1):ctx.box.N(2);

% patternsearch options
options = optimoptions("patternsearch", ...
    "Display", "off", ...
    "UseCompletePoll", true, ...
    "MeshTolerance", 1e-5);

% Parallel pool
pool = gcp();
wait(parfevalOnAll(pool, @addpath, 0, genpath(repoRoot)));  % let the workers see the functions

%% Enumerate N in parallel, solve the 2D (x1, h) problem for each
[X, Jbest] = run_shape_search(N_list, lb, ub, options, params, ctx, n_starts, opts.quiet);

%% The circular port, as a candidate on equal terms
% The star family cannot degenerate into a circle, because C2 needs a tip at
% least a few cells tall; the optimizer pressing against that floor is the
% search asking for a cylinder. So the cylinder is optimized separately, over
% its one variable, and ranked next to the stars. It is neither forced nor
% penalized: if it wins, it wins.
[x1_cyl, J_cyl] = run_cylinder_search(params, ctx, opts.quiet);

%% Ranking over all N, RE-OPTIMIZED on the fine grid

params_fine = params;
params_fine.mdf.grid_divisions = grid_fine;
ctx_fine = ctx;
ctx_fine.h_min = ctx.K.C2.lo / grid_fine;

lb_fine = [ctx_fine.box.x1(1); ctx_fine.h_min];
ub_fine = [ctx_fine.box.x1(2); 1 - ctx_fine.box.x1(1)];

rank_table = build_ranking(N_list, X, Jbest, x1_cyl, J_cyl, ...
    params_fine, ctx_fine, lb_fine, ub_fine, options, opts.quiet);

if isempty(rank_table)
    shapes = struct([]);
    if ~opts.quiet
        fprintf("Phase A found no feasible shape for %s.\n", ctx.oxidizer);
    end
    return
end

rank_table = rank_table(1:min(n_report, height(rank_table)), :);

if ~opts.quiet
    n_show = height(rank_table);
    fprintf("\n RANKING - top %d shapes over all N (%s, fine grid %d)\n", ...
        n_show, ctx.oxidizer, grid_fine);
    fprintf("%s\n", repmat('-', 1, 100));
    fprintf(" %-4s %-9s %-3s %7s %8s %10s %7s %8s %7s %8s %7s\n", ...
        "rank", "geometry", "N", "x1", "h", "Isp_load", "OF_med", "sigma", ...
        "drift", "mdot_ox", "Gox0");
    fprintf("%s\n", repmat('-', 1, 100));
    for i = 1:n_show
        r = rank_table(i, :);
        fprintf(" %-4d %-9s %-3d %7.4f %8.5f %10.2f %7.3f %8.4f %7.3f %8.2f %7.0f\n", ...
            i, r.geometry, r.N, r.x1, r.h, r.Isp_load, r.OF_med, r.sigma, ...
            r.drift, r.mdot_ox, r.Gox0);
    end
    fprintf("%s\n", repmat('-', 1, 100));

    best = rank_table(1, :);
    fprintf("\nBest shape (%s): %s\n", ctx.oxidizer, best.geometry);
    if best.geometry == "star"
        fprintf("\tN            = %d\n", best.N);
        fprintf("\th  (tip)     = %.5f R_c   (C2 floor: %.5f)\n", best.h, ctx_fine.h_min);
        fprintf("\tx2 (re/R_c)  = %.4f\n", best.x2);
    end
    fprintf("\tx1 (port)    = %.4f\n", best.x1);
    fprintf("\tIsp_load     = %.2f s   <- objective, from I_tot/(g0 m_load)\n", best.Isp_load);
    fprintf("\tIsp_med      = %.2f s   (sliver factor %.4f)\n", ...
        best.Isp_med, (best.OF_med + 1)/(best.OF_med + 1/(1 - best.sigma)));
    fprintf("\tO/F mean     = %.3f    (%.3f -> %.3f)\n", ...
        best.OF_med, best.of0, best.of_end);
    fprintf("\tdrift        = %.3f  useful burn / %.3f  to burnout\n", ...
        best.drift, best.drift_full);
    fprintf("\tsliver sigma = %.4f\n", best.sigma);
    fprintf("\tmdot_ox      = %.3f kg/s  <- computed from C7, not searched\n", best.mdot_ox);
    fprintf("\tG_ox(0)      = %.0f kg/m2 s  (C4: [%.0f, %.0f], C10: [%.0f, %.0f])\n", ...
        best.Gox0, ctx.K.C4.lo, ctx.K.C4.hi, ctx.K.C10.lo, ctx.K.C10.hi);
    fprintf("\tGtilde       = %.3f\n", best.Gtilde);
    fprintf("\tpredicted    : 2R_c = %.3f m, L = %.3f m, m_load = %.0f kg\n", ...
        2*best.R_c, best.L, best.m_load);
end

best = rank_table(1, :);

%% Save the ranking
% The whole ranking is written, with a run stamp and the constraint
% configuration in force, so phase B can read the BEST row of the BEST run
if opts.write_csv
    run_id = string(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss"));
    out = rank_table;
    out.run_id = repmat(run_id, height(out), 1);
    out.rank = (1:height(out)).';
    out.grid_divisions = repmat(grid_fine, height(out), 1);
    out.oxidizer = repmat(string(ctx.oxidizer), height(out), 1);
    out.Gox_lo = repmat(ctx.K.C10.lo, height(out), 1);
    out.Gox_hi = repmat(ctx.K.C10.hi, height(out), 1);
    out.of_lo = repmat(ctx.K.C5.lo, height(out), 1);
    out.of_hi = repmat(ctx.K.C5.hi, height(out), 1);
    out.thrust_mode = repmat(string(ctx.K.C7.value), height(out), 1);

    archive_incompatible_csv(resultFile, out, repoRoot);
    if isfile(resultFile)
        writetable(out, resultFile, "WriteMode", "append");
    else
        writetable(out, resultFile);
    end
    if ~opts.quiet
        fprintf("\nRanking written to %s (run %s)\n", resultFile, run_id);
    end
end

%% The ranking, as a struct array
shapes = table2struct(rank_table).';
for i = 1:numel(shapes)
    shapes(i).oxidizer = string(ctx.oxidizer);
    shapes(i).rank = i;
end

if ~opts.do_plots
    return
end

%% Plots
figure("Name", "Phase A - best shape", "Color", "w");
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");

R_c = params_fine.engine.ext_diameter / 2;
th = linspace(0, 2*pi, 400);
if best.geometry == "cylinder"
    lookup = cylinderLookup(best.x1, params_fine);
    r_port = best.x1 * R_c;
    port = [r_port*cos(th).', r_port*sin(th).'];
else
    meshdata.inner_diameter = 2 * best.x1 * R_c;
    meshdata.outer_diameter = 2 * best.x2 * R_c;
    meshdata.n_tips = best.N;
    lookup = build_shape_lookup(meshdata, params_fine);
    half = make_mesh0("star", meshdata, 800, "cartesian");
    port = [half; flipud([half(:,1), -half(:,2)])];
end
Phi = (lookup.Ap / R_c^2).^params.fuel.n_rf ./ max(lookup.perim / R_c, realmin);

% Grain cross-section
nexttile
fill(R_c*cos(th), R_c*sin(th), [0.82 0.72 0.47], "EdgeColor", "k", "LineWidth", 1.3);
hold on
fill(port(:,1), port(:,2), "w", "EdgeColor", [0.20 0.20 0.20], "LineWidth", 1.1);
axis equal
axis(R_c * [-1.1 1.1 -1.1 1.1])
set(gca, "XTick", [], "YTick", [])
if best.geometry == "cylinder"
    title(sprintf("circular port, x1 = %.3f", best.x1));
else
    title(sprintf("star N = %d, x1 = %.3f, h = %.4f", best.N, best.x1, best.h));
end

% Shape function, which is the O/F history up to the level lambda
nexttile
plot(lookup.b / R_c, Phi, "LineWidth", 1.6, "Color", [0 0.45 0.74]);
grid on
xlabel("$\tilde b$ [-]", "Interpreter", "latex");
ylabel("$\tilde \Phi = \tilde A_p^n / \tilde P_b$ [-]", "Interpreter", "latex");
title(sprintf("O/F drift %.3f, sliver %.4f", best.drift, best.sigma));

% Objective against tip count
figure("Name", "Phase A - objective vs N", "Color", "w");
plot(N_list, -Jbest, "o-", "LineWidth", 1.5);
hold on
plot(best.N, best.Isp_load, "p", "MarkerSize", 14, ...
    "MarkerFaceColor", [0.85 0.33 0.10], "MarkerEdgeColor", "k");
grid on
xlabel("number of tips N [-]");
ylabel("best I_{sp} on loaded mass [s]");
title("Phase A objective vs tip count");

end

%% FUNCTIONS

function opts = default_opt(opts, name, value)
    % default_opt
    % Fill in one optional field of the options struct.
    % INPUT
    %   opts  / struct / 1x1   Options
    %   name  / string / 1x1   Field name
    %   value / any    / 1x1   Default value
    % OUTPUT
    %   opts  / struct / 1x1   Options with the field guaranteed present

    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = value;
    end
end

function [X, Jbest] = run_shape_search(N_list, lb, ub, options, params, ctx, n_starts, quiet)
    % run_shape_search
    % Solve the fixed-N two-dimensional shape problem for every N in parallel,
    % from several starting points, and report live progress from the client
    % through a DataQueue.
    % INPUT
    %   N_list   / double / 1xM   Tip counts to enumerate [-]
    %   lb, ub   / double / 2x1   Bounds on [x1; h] [-]
    %   options  / optim  / 1x1   patternsearch options
    %   params   / struct / 1x1   Combustion parameters
    %   ctx      / struct / 1x1   Evaluation context (engineContext)
    %   n_starts / double / 1x1   Starting points per tip count [-]
    %   quiet    / logical/ 1x1   Suppress the progress report
    % OUTPUT
    %   X        / double / Mx2   Best (x1, h) per N
    %   Jbest    / double / Mx1   Best cost per N [s]

    n_N = numel(N_list);
    X = nan(n_N, 2);
    Jbest = inf(n_N, 1);

    % Progress indicator: afterEach runs on the client each time a worker sends
    q = parallel.pool.DataQueue;
    done = 0;
    afterEach(q, @report_progress);

    if ~quiet
        fprintf("Phase A (%s): %d tip counts in parallel, %d starts each...\n", ...
            ctx.oxidizer, n_N, n_starts);
    end

    parfor k = 1:n_N
        N = N_list(k);

        % Multi-start
        x0_list = starting_points(N, lb, ub, n_starts);

        fun = @(z) shapeCostFunction(z(1), z(2), N, params, ctx);

        best_x = [NaN, NaN];
        best_J = Inf;
        for s = 1:size(x0_list, 1)
            % Containment C3 as the only linear inequality left:
            % x1/cos(pi/N) + h <= 1
            Aineq = [1/cos(pi/N), 1];
            [x, fval] = patternsearch(fun, x0_list(s,:), Aineq, 1, [], [], ...
                lb, ub, [], options);
            if fval < best_J
                best_J = fval;
                best_x = x;
            end
        end

        X(k,:) = best_x;
        Jbest(k) = best_J;
        send(q, k);
    end

    function report_progress(~)
        done = done + 1;
        if ~quiet
            fprintf("  ... %2d/%d done (%.0f%%)\n", done, n_N, 100*done/n_N);
        end
    end
end

function x0_list = starting_points(N, lb, ub, n_starts)
    % starting_points
    % Starting points for one tip count, spread over the port fraction: a bland
    % tip at a small port, a bland tip at a mid port, and a deep star.
    % INPUT
    %   N        / double / 1x1   Number of tips [-]
    %   lb, ub   / double / 2x1   Bounds on [x1; h] [-]
    %   n_starts / double / 1x1   How many starting points [-]
    % OUTPUT
    %   x0_list  / double / Sx2   Starting points, rows of [x1, h]

    c = cos(pi/N);
    candidates = [ ...
        0.22,             max(lb(2), 0.005); ...   % bland tip, small port
        0.45,             max(lb(2), 0.020); ...   % bland tip, mid port
        0.5*c,            0.9 - 0.5];              % deep star (the old start)

    n_starts = min(n_starts, size(candidates, 1));
    x0_list = candidates(1:n_starts, :);

    % Clip into the box and onto the containment face
    x0_list(:,1) = min(max(x0_list(:,1), lb(1)), ub(1));
    x0_list(:,2) = min(max(x0_list(:,2), lb(2)), ub(2));
    over = x0_list(:,1)/c + x0_list(:,2) > 1;
    x0_list(over, 2) = max(lb(2), 1 - x0_list(over,1)/c - 1e-6);
end

function archive_incompatible_csv(resultFile, out, repoRoot)
    % archive_incompatible_csv
    % Move an existing result file out of the way when its columns no longer
    % match what this run writes. Old runs are archived, never deleted.
    % INPUT
    %   resultFile / string / 1x1   Path of the CSV
    %   out        / table  / Nx?   Table about to be written
    %   repoRoot   / string / 1x1   Repository root
    % OUTPUT
    %   None (may rename the existing file into Dumpster/)

    if ~isfile(resultFile)
        return
    end
    old = readtable(resultFile);
    if isequal(string(old.Properties.VariableNames), string(out.Properties.VariableNames))
        return
    end

    dumpster = fullfile(repoRoot, "Dumpster");
    if ~isfolder(dumpster)
        mkdir(dumpster);
    end
    [~, name, ext] = fileparts(resultFile);
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    archived = fullfile(dumpster, name + "_legacy_" + stamp + ext);
    movefile(resultFile, archived);
    fprintf("Existing %s%s had the old columns: archived to %s\n", name, ext, archived);
end

function [x1_opt, J_opt] = run_cylinder_search(params, ctx, quiet)
    % run_cylinder_search
    % Optimize the circular port over its single variable, the port radius
    % ratio. Closed-form geometry, so this is a one-dimensional problem that
    % costs a fraction of a single star evaluation.
    % INPUT
    %   params / struct / 1x1   Combustion parameters
    %   ctx    / struct / 1x1   Evaluation context (engineContext)
    %   quiet  / logical/ 1x1   Suppress the report
    % OUTPUT
    %   x1_opt / double / 1x1   Best port radius ratio [-]
    %   J_opt  / double / 1x1   Best cost [s]

    lo = ctx.box.x1(1);
    hi = ctx.box.x1(2);

    % Coarse scan first: C4 carves an admissible interval out of the box, and
    % fminbnd on its own would happily converge inside the infeasible region.
    % The scan runs on the same coarse MDF grid as the star search, and in
    % parallel, because each point is now a full MDF solve.
    n_scan = 40;
    x_scan = linspace(lo, hi, n_scan);
    J_scan = inf(1, n_scan);
    parfor i = 1:n_scan
        J_scan(i) = cylinderCostFunction(x_scan(i), params, ctx);
    end

    [J_opt, ibest] = min(J_scan);
    if ~isfinite(J_opt)
        x1_opt = NaN;
        if ~quiet
            fprintf("Cylinder: no admissible port radius.\n");
        end
        return
    end
    x1_opt = x_scan(ibest);

    % Refine inside the bracket around the best scan point
    a = x_scan(max(1, ibest - 1));
    c = x_scan(min(n_scan, ibest + 1));
    if c > a
        [x_ref, J_ref] = fminbnd(@(x) cylinderCostFunction(x, params, ctx), a, c, ...
            optimset("Display", "off", "TolX", 1e-4));
        if isfinite(J_ref) && J_ref < J_opt
            x1_opt = x_ref;
            J_opt = J_ref;
        end
    end

    if ~quiet
        fprintf("Cylinder: x1 = %.4f, Isp_load = %.2f s\n", x1_opt, -J_opt);
    end
end

function T = build_ranking(N_list, X, Jbest, x1_cyl, J_cyl, params_fine, ...
        ctx_fine, lb, ub, options, quiet)
    % build_ranking
    % RE-OPTIMIZE every per-N optimum on the fine grid, re-optimize the circular
    % port too, and sort the lot by merit.
    %
    % Re-optimizing rather than re-scoring is not a refinement, it is a
    % correctness requirement: C2 floors h at 5/grid_divisions, the search sits
    % exactly on that floor, and the floor halves when the grid doubles. A
    % re-scored star would be ranked while stuck at twice the floor the ranking
    % declares, and since the cylinder carries no grid-derived bound the error
    % would fall on one family only.
    %
    % A shape that stops being admissible on the fine grid is REPORTED, not
    % dropped in silence: it disappearing from the ranking without a word is
    % exactly the kind of thing that hides a constraint problem.
    % INPUT
    %   N_list      / double / 1xM   Tip counts [-]
    %   X           / double / Mx2   Best (x1, h) per N on the coarse grid [-]
    %   Jbest       / double / Mx1   Coarse-grid cost per N [s]
    %   x1_cyl      / double / 1x1   Best circular port radius ratio [-]
    %   J_cyl       / double / 1x1   Its cost [s]
    %   params_fine / struct / 1x1   Parameters on the fine grid
    %   ctx_fine    / struct / 1x1   Context with the fine-grid C2 floor
    %   lb, ub      / double / 2x1   Bounds with the FINE-grid h floor [-]
    %   options     / optim  / 1x1   patternsearch options
    %   quiet       / logical/ 1x1   Suppress the report
    % OUTPUT
    %   T           / table  / Kx?   Feasible shapes, best first

    n_N = numel(N_list);
    rows = cell(n_N + 1, 1);
    moved = nan(n_N, 4);        % [h_before, h_after, J_before, J_after]
    dropped = strings(1, 0);

    % A short leash: the restart is already close, so this is a descent onto the
    % new floor, not a fresh search
    opt_fine = options;
    opt_fine.MaxIterations = 40;

    for k = 1:n_N
        if ~isfinite(Jbest(k)) || any(isnan(X(k,:)))
            continue
        end
        N = N_list(k);

        % Clip the coarse optimum into the fine-grid box before restarting
        x0 = [min(max(X(k,1), lb(1)), ub(1)), min(max(X(k,2), lb(2)), ub(2))];
        Aineq = [1/cos(pi/N), 1];          % C3 containment

        fun = @(z) shapeCostFunction(z(1), z(2), N, params_fine, ctx_fine);
        [x, ~] = patternsearch(fun, x0, Aineq, 1, [], [], lb, ub, [], opt_fine);

        [Jf, info] = fun(x);
        moved(k,:) = [X(k,2), x(2), -Jbest(k), -Jf];

        if ~isfinite(Jf) || ~info.feasible
            dropped(end+1) = sprintf("N = %d (%s)", N, info.fail); %#ok<AGROW>
            continue
        end
        rows{k} = ranking_row("star", N, x(1), x(2), Jf, info);
    end

    % The cylinder is re-optimized on the same fine grid, for the same reason
    if isfinite(J_cyl) && ~isnan(x1_cyl)
        lo = max(ctx_fine.box.x1(1), x1_cyl*0.85);
        hi = min(ctx_fine.box.x1(2), x1_cyl*1.15);
        x1_ref = x1_cyl;
        if hi > lo
            x1_ref = fminbnd(@(x) cylinderCostFunction(x, params_fine, ctx_fine), ...
                lo, hi, optimset("Display", "off", "TolX", 1e-5));
        end
        [Jf, info] = cylinderCostFunction(x1_ref, params_fine, ctx_fine);
        [Jf0, info0] = cylinderCostFunction(x1_cyl, params_fine, ctx_fine);
        if ~isfinite(Jf) || ~info.feasible || (isfinite(Jf0) && Jf0 < Jf)
            Jf = Jf0;
            info = info0;
            x1_ref = x1_cyl;
        end
        if isfinite(Jf) && info.feasible
            rows{end} = ranking_row("cylinder", 0, x1_ref, 0, Jf, info);
        else
            dropped(end+1) = sprintf("cylinder (%s)", info.fail);
        end
    end

    % Report what the refinement actually did, so a re-optimization that is not
    % running can be spotted: h would simply not have moved
    if ~quiet
        ok_moved = ~isnan(moved(:,1));
        if any(ok_moved)
            fprintf("\n Fine-grid re-optimization (C2 floor on h is now %.5f):\n", ...
                ctx_fine.h_min);
            fprintf("  %-4s %10s %10s %10s %10s\n", ...
                "N", "h before", "h after", "Isp before", "Isp after");
            for k = find(ok_moved).'
                fprintf("  %-4d %10.5f %10.5f %10.2f %10.2f\n", ...
                    N_list(k), moved(k,1), moved(k,2), moved(k,3), moved(k,4));
            end
        end
        if ~isempty(dropped)
            fprintf(" Shapes dropped on the fine grid: %s\n", strjoin(dropped, ", "));
        end
    end

    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows)
        T = table();
        return
    end
    T = vertcat(rows{:});
    T = sortrows(T, "Isp_load", "descend");
end

function row = ranking_row(geometry, N, x1, h, Jf, info)
    % ranking_row
    % One row of the phase A ranking, identical in shape for both families.
    % INPUT
    %   geometry / string / 1x1   "star" or "cylinder"
    %   N        / double / 1x1   Tip count, 0 for the cylinder [-]
    %   x1, h    / double / 1x1   Port ratio and tip height [-]
    %   Jf       / double / 1x1   Cost [s]
    %   info     / struct / 1x1   Diagnostics from the cost function
    % OUTPUT
    %   row      / table  / 1x?   Ranking row

    row = table(string(geometry), N, x1, h, info.x2, -Jf, info.Isp_med, ...
        info.OF_med, info.sigma, info.drift, info.drift_full, info.Gtilde, ...
        info.lambda, info.Phi0, info.of0, info.of_end, info.b_end, ...
        info.mdot_ox, info.Gox0, info.R_c, info.L, info.m_load, info.I_tot, ...
        'VariableNames', {'geometry','N','x1','h','x2','Isp_load','Isp_med', ...
        'OF_med','sigma','drift','drift_full','Gtilde','lambda','Phi0', ...
        'of0','of_end','b_end','mdot_ox','Gox0','R_c','L','m_load','I_tot'});
end
