function design = phaseB(shapes, thermo, params, C, opts)
% phaseB
% PHASE B - engine sizing on top of the phase A shape ranking.
%
% The shape and the oxidizer flow both come from phase A, which consumed the
% single degree of freedom of the problem when it picked the O/F level lambda
% (four unknowns, three equations). Phase B is therefore not a search any more:
% it computes the exact design with the ODE model and refines the flow locally.
% For each candidate flow, sizeEngine solves C6, C7, C8 for R_c, L and A_t, and
% the design is scored by the specific impulse on the LOADED mass,
%   J = -I_tot/(g0*m_load),   m_load = mdot_ox*t_b + rho_f*L*(pi*R_c^2 - Ap(0))
% which is the mass that leaves the ground, sliver included. There are no
% penalties: C5, C9, C10, C11, C12 and C13 are boolean gates.
%
% The shapes of the phase A ranking are tried in order until one produces a
% design that passes every gate. Every candidate is evaluated with the SAME
% exact ODE model as the final validation run.
% INPUT
%   shapes / struct / 1xK   Ranked shapes from phaseA
%   thermo / struct / 1x1   Oxidizer thermochemistry, from get_thermo
%   params / struct / 1x1   Combustion parameters (from combustion_params())
%   C      / struct / 1xM   Constraint table from optimizationConstraints
%   opts   / struct / 1x1   Optional: do_plots (false), quiet (false),
%                           n_shapes (5), n_sweep (7), bracket (0.30),
%                           validate (false)
% OUTPUT
%   design / struct / 1x1   ok, failTag, and on success the sized engine S,
%                           the shape used and the sweep record

if nargin < 5 || isempty(opts)
    opts = struct();
end
opts = default_opt(opts, "do_plots", false);
opts = default_opt(opts, "quiet", false);
opts = default_opt(opts, "n_shapes", 5);       % [-] shapes to try before giving up
opts = default_opt(opts, "n_sweep", 7);        % [-] flow grid points in the bracket
opts = default_opt(opts, "validate", false);   % re-run the winner from scratch

design = struct("ok", false, "failTag", "", "oxidizer", string(thermo.name));

repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));

params.geometry.type = "star";
params.mdf.grid_divisions = 700;     % fine grid, built once per shape and rescaled
params.mdf.perimeter_from_area = false;

% Shared context: constraints, thermochemistry, Isp and c_eff curves
ctx = engineContext(thermo, params, C);

if isempty(shapes)
    design.failTag = "phase A returned no shape";
    return
end
n_try = min(opts.n_shapes, numel(shapes));
if ~opts.quiet
    fprintf("\nPhase B (%s): %d candidate shapes\n", ctx.oxidizer, n_try);
    for i = 1:n_try
        fprintf("  %d) %-8s N = %2d, x1 = %.4f, h = %.5f, Isp_load = %.2f s\n", ...
            i, shapes(i).geometry, shapes(i).N, shapes(i).x1, shapes(i).h, ...
            shapes(i).Isp_load);
    end
end

%% Sizing configuration
cfg = struct();
cfg.n = params.fuel.n_rf;
cfg.a = params.fuel.a_rf;
cfg.rho_f = params.fuel.rho_f;
cfg.web_min = ctx.K.C12.lo;                   % [m]   C12
cfg.p_min = params.combustion.p_min;          % [Pa]  lookup extent
cfg.p_max = params.combustion.p_max;          % [Pa]  lookup extent
% The horizon must exceed the burn-time target, or ode113 would cap burn_time
% at t_target and the C6 equation would be one-way
cfg.time_output = 0:1.0:1.5*ctx.K.C6.lo;      % [s]
cfg.fine_ode = params.time.fine_ode;

% Oxidizer-flow bracket. Phase A already CONSUMED the single degree of freedom
% of the problem when it picked the O/F level (four unknowns, three equations),
% so phase B is not a search: it is an exact calculation plus a local refinement
% around the flow phase A computed. A wide bracket would send the optimizer back
% to exploring flows the thrust requirement does not allow, which is exactly the
% mistake the old existential form of C4 made.
half_width = ctx.box.mdot_bracket;             % [-] +/-30 % by default

% Parallel pool
pool = gcp();
wait(parfevalOnAll(pool, @addpath, 0, genpath(repoRoot)));

%% Try the phase A shapes in order until one passes every gate
winner = [];
fail_tags = strings(1, n_try);
for s = 1:n_try
    if ~opts.quiet
        fprintf("\n===== Shape %d/%d: %s, N = %d, x1 = %.4f, h = %.5f =====\n", ...
            s, n_try, shapes(s).geometry, shapes(s).N, shapes(s).x1, shapes(s).h);
    end

    % Normalized geometry at R_c = 1, built once and reused by rescaling.
    % The circular port has a closed form; the star needs the MDF field.
    params_norm = params;
    params_norm.engine.ext_diameter = 2.0;         % R_c = 1 m
    if shapes(s).geometry == "cylinder"
        params_norm.geometry.type = "cylinder";
        lookupN = cylinderLookup(shapes(s).x1, params_norm);
    else
        meshdata = struct();
        meshdata.inner_diameter = 2 * shapes(s).x1;
        meshdata.outer_diameter = 2 * shapes(s).x2;
        meshdata.n_tips = shapes(s).N;
        lookupN = build_shape_lookup(meshdata, params_norm);
    end

    cfg_s = cfg;
    cfg_s.x2 = shapes(s).x2;                       % [-] for the C12 burn stop

    % Bracket around the flow phase A computed for THIS shape
    m0 = shapes(s).mdot_ox;
    mdot_grid = linspace(m0*(1 - half_width), m0*(1 + half_width), opts.n_sweep).';
    if ~opts.quiet
        fprintf("  mdot_ox from phase A: %.3f kg/s, bracket [%.2f, %.2f]\n", ...
            m0, mdot_grid(1), mdot_grid(end));
    end

    [best_s, sweep_s] = optimize_mdot(mdot_grid, ctx, lookupN, cfg_s, opts.quiet);
    if ~opts.quiet
        print_sweep(mdot_grid, sweep_s);
    end

    if ~isempty(best_s)
        winner = best_s;
        winner.shape = shapes(s);
        winner.lookupN = lookupN;
        winner.cfg = cfg_s;
        winner.sweep = sweep_s;
        winner.mdot_grid = mdot_grid;
        if ~opts.quiet
            fprintf("\nShape %d passes every gate.\n", s);
        end
        break
    end
    fail_tags(s) = collect_tags(sweep_s);
    if ~opts.quiet
        fprintf("\nShape %d has no admissible oxidizer flow (%s), trying the next one.\n", ...
            s, fail_tags(s));
    end
end

if isempty(winner)
    % Never discard an oxidizer in silence: report which gates refused it
    design.ok = false;
    design.failTag = "no admissible design; gates hit: " + ...
        strjoin(unique(fail_tags(strlength(fail_tags) > 0)), ", ");
    if ~opts.quiet
        fprintf("\nPhase B (%s) failed: %s\n", ctx.oxidizer, design.failTag);
    end
    return
end

S = winner.S;
shape = winner.shape;

design.ok = true;
design.failTag = "OK";
design.S = S;
design.shape = shape;
design.sweep = winner.sweep;
design.mdot_grid = winner.mdot_grid;
design.ctx = ctx;
% Enough to re-run the exact model on this design from the outside, which is
% what verifySolver needs for the three sweeps of section 7.1
design.lookupN = winner.lookupN;
design.cfg = winner.cfg;


%% Present the design: characteristics, figures, independent re-run.
% The three switches are independent: quiet controls the printing, do_plots
% the figures, validate the rebuild through run_combustion_case. None of them
% gates the others, and reportDesign can be called later on a design that has
% already been sized, without paying for the sizing again.
design = reportDesign(design, params, struct( ...
    "quiet", opts.quiet, "do_plots", opts.do_plots, "validate", opts.validate));

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

function tag = collect_tags(sweep)
    % collect_tags
    % Gather the distinct gate tags of a whole sweep, so a shape that fails is
    % reported with the reason and not just dropped.
    % INPUT
    %   sweep / struct / 1xM   Records from optimize_mdot
    % OUTPUT
    %   tag   / string / 1x1   Distinct tags, comma separated

    tags = strings(1, numel(sweep));
    for j = 1:numel(sweep)
        tags(j) = string(sweep(j).tag);
    end
    tags = unique(tags(strlength(tags) > 0));
    tag = strjoin(tags, ", ");
end

function shapes = load_shape_ranking(shapeFile, n_shapes)
    % load_shape_ranking
    % Read the phase A ranking from disk and return the best shapes of the most
    % recent run, best first. Kept for running phase B on its own, without
    % main_optimization: reading the last row appended would mean sizing
    % whatever run happened to finish last, not the best shape found.
    % INPUT
    %   shapeFile / string / 1x1   Path of the phase A CSV
    %   n_shapes  / double / 1x1   How many shapes to keep [-]
    % OUTPUT
    %   shapes    / table  / Kx?   Candidate shapes, best first

    if ~isfile(shapeFile)
        error("Shape ranking not found at %s. Run shapeOptimization first.", shapeFile);
    end
    T = readtable(shapeFile, "TextType", "string");
    if ~ismember("run_id", T.Properties.VariableNames)
        error(['%s has the old columns. Re-run shapeOptimization: phase B needs ' ...
            'the ranking format.'], shapeFile);
    end

    last_run = T.run_id(find(T.run_id == max(T.run_id), 1, "last"));
    T = T(T.run_id == last_run, :);
    T = sortrows(T, "Isp_load", "descend");
    shapes = T(1:min(n_shapes, height(T)), :);
end

function [best, sweep] = optimize_mdot(mdot_grid, ctx, lookupN, cfg, quiet)
    % optimize_mdot
    % Sweep the oxidizer flow in parallel, size the engine at each point, apply
    % the gates, then refine around the best admissible point.
    % INPUT
    %   mdot_grid / double / Mx1   Oxidizer flows to try [kg/s]
    %   ctx       / struct / 1x1   Context from engineContext
    %   lookupN   / struct / 1x1   Normalized geometry lookup at R_c = 1
    %   cfg       / struct / 1x1   Sizing configuration for sizeEngine
    %   quiet     / logical/ 1x1   Suppress the progress report
    % OUTPUT
    %   best      / struct / 1x1   Winning design, empty if no point passes
    %   sweep     / struct / 1xM   One record per grid point

    n_m = numel(mdot_grid);
    sweep = repmat(struct("ok", false, "feasible", false, "Isp_load", NaN, ...
        "tag", "", "S", []), 1, n_m);

    % Progress indicator: afterEach runs on the client each time a worker sends
    q = parallel.pool.DataQueue;
    done = 0;
    afterEach(q, @report_progress);
    if ~quiet
        fprintf("Sizing %d oxidizer flows in parallel...\n", n_m);
    end

    parfor j = 1:n_m
        S = sizeEngine(mdot_grid(j), ctx, lookupN, cfg);
        rec = struct("ok", S.ok, "feasible", false, "Isp_load", NaN, ...
            "tag", "", "S", []);
        if S.ok
            [feasible, tag] = feas_check(S, ctx);
            rec.feasible = feasible;
            rec.tag = tag;
            rec.Isp_load = S.Isp_load;
            rec.S = S;
        else
            rec.tag = "sizing: " + S.err;
        end
        sweep(j) = rec;
        send(q, j);
    end

    % Best admissible grid point
    score = -inf(1, n_m);
    for j = 1:n_m
        if sweep(j).ok && sweep(j).feasible
            score(j) = sweep(j).Isp_load;
        end
    end
    [best_score, jbest] = max(score);
    if ~isfinite(best_score)
        best = [];
        return
    end

    % Refine inside the admissible bracket around it
    lo = mdot_grid(max(1, jbest - 1));
    hi = mdot_grid(min(n_m, jbest + 1));
    best = struct("S", sweep(jbest).S, "tag", sweep(jbest).tag);
    if hi > lo
        obj = @(m) neg_Isp_load(m, ctx, lookupN, cfg);
        m_opt = fminbnd(obj, lo, hi, optimset("Display", "off", "TolX", 0.02));
        S_opt = sizeEngine(m_opt, ctx, lookupN, cfg);
        if S_opt.ok
            [feasible, tag] = feas_check(S_opt, ctx);
            % Only accept the refinement if it is admissible AND better: near
            % the edge of a gate the refined point can fall just outside
            if feasible && S_opt.Isp_load > best.S.Isp_load
                best.S = S_opt;
                best.tag = tag;
            end
        end
    end

    function report_progress(~)
        done = done + 1;
        if ~quiet
            fprintf("  ... %2d/%d done (%.0f%%)\n", done, n_m, 100*done/n_m);
        end
    end
end

function J = neg_Isp_load(mdot_ox, ctx, lookupN, cfg)
    % neg_Isp_load
    % Objective of phase B: minimize minus the specific impulse on the LOADED
    % mass. No penalties, no weights; inadmissible points are simply rejected.
    % INPUT
    %   mdot_ox / double / 1x1   Oxidizer mass flow [kg/s]
    %   ctx     / struct / 1x1   Context from engineContext
    %   lookupN / struct / 1x1   Normalized geometry lookup
    %   cfg     / struct / 1x1   Sizing configuration
    % OUTPUT
    %   J       / double / 1x1   -Isp_load [s], +Inf when the design fails

    S = sizeEngine(mdot_ox, ctx, lookupN, cfg);
    if ~S.ok || ~isfinite(S.Isp_load)
        J = Inf;
        return
    end
    [feasible, ~] = feas_check(S, ctx);
    if ~feasible
        J = Inf;
        return
    end
    J = -S.Isp_load;
end

function [feasible, tag] = feas_check(S, ctx)
    % feas_check
    % Apply the phase B gates, all read from the constraint table: no threshold
    % is written here. Every one is boolean, none is weighted.
    % INPUT
    %   S   / struct / 1x1   Sized engine from sizeEngine
    %   ctx / struct / 1x1   Context carrying the constraint table
    % OUTPUT
    %   feasible / logical / 1x1   True when every gate passes
    %   tag      / string  / 1x1   "OK", or the ids that failed

    K = ctx.K;
    if ~S.ok
        feasible = false;
        tag = "sizing";
        return
    end

    tol_of = 1e-6;    % [-] the O/F event terminates exactly on the bound
    ok_C5 = S.of_min_hist >= K.C5.lo - tol_of && S.of_max_hist <= K.C5.hi + tol_of;
    ok_C9 = S.thrust0 >= K.C9.lo && S.thrust0 <= K.C9.hi;
    ok_C10 = S.Gox0 >= K.C10.lo && S.Gox0 <= K.C10.hi;
    ok_C11 = S.Gox_min >= K.C11.lo;
    ok_C12 = S.web_residual >= K.C12.lo - 1e-9;
    ok_C13 = 2*S.R_c <= K.C13.hi(1) && S.L <= K.C13.hi(2);

    feasible = ok_C5 && ok_C9 && ok_C10 && ok_C11 && ok_C12 && ok_C13;

    tag = "";
    if ~ok_C5,  tag = tag + "C5 ";  end
    if ~ok_C9,  tag = tag + "C9 ";  end
    if ~ok_C10, tag = tag + "C10 "; end
    if ~ok_C11, tag = tag + "C11 "; end
    if ~ok_C12, tag = tag + "C12 "; end
    if ~ok_C13, tag = tag + "C13 "; end
    if feasible
        tag = "OK";
    else
        tag = strtrim(tag);
    end
end

function print_sweep(mdot_grid, sweep)
    % print_sweep
    % Print the oxidizer-flow sweep, with the sized engine, the objective and
    % the failing gate of every point.
    % INPUT
    %   mdot_grid / double / Mx1   Oxidizer flows [kg/s]
    %   sweep     / struct / 1xM   Records from optimize_mdot
    % OUTPUT
    %   None (prints to stdout)

    fprintf("\n %6s %9s %8s %7s %7s %8s %8s %8s  %s\n", ...
        "mdot", "Isp_load", "F0[kN]", "2R_c[m]", "L[m]", "Gox0", "t_b[s]", "sigma", "gate");
    fprintf("%s\n", repmat('-', 1, 88));
    for j = 1:numel(mdot_grid)
        rec = sweep(j);
        if ~rec.ok
            fprintf(" %6.2f %9s %8s %7s %7s %8s %8s %8s  %s\n", ...
                mdot_grid(j), "-", "-", "-", "-", "-", "-", "-", rec.tag);
            continue
        end
        S = rec.S;
        fprintf(" %6.2f %9.2f %8.1f %7.3f %7.3f %8.0f %8.1f %8.4f  %s\n", ...
            mdot_grid(j), S.Isp_load, S.thrust0*1e-3, 2*S.R_c, S.L, ...
            S.Gox0, S.burn_time, S.sigma, rec.tag);
    end
    fprintf("%s\n", repmat('-', 1, 88));
end

function compare_row(name, a, b)
    % compare_row
    % One line of the sizing-model against full-rebuild comparison.
    % INPUT
    %   name / string / 1x1   Quantity name
    %   a    / double / 1x1   Value from the sizing model
    %   b    / double / 1x1   Value from the independent rebuild
    % OUTPUT
    %   None (prints to stdout)

    fprintf(" %-18s %12.4f %12.4f %7.2f%%\n", name, a, b, 100*(b/a - 1));
end
