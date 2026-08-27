% MAIN_OPTIMIZATION - single entry point of the two-phase engine optimization.
%
% The oxidizer parametrizes the WHOLE A -> B chain: every
% pair has its own peak of Isp(O/F), so it moves the best O/F level, hence
% mdot_ox (both through c_eff and through OF/(1+OF)), hence G_ox(0), hence the
% admissible window on the port fraction, hence potentially the optimal shape
% itself. Re-running phase B on shapes chosen for a different oxidizer would
% therefore answer the wrong question, and the loop below is over the full
% chain.
%
% The loop is SERIAL on purpose: the parfor lives inside phaseA, on the tip
% count, where sixteen short iterations balance far better than four or five
% long ones. Nesting the two would only fight over the same workers.
%
% An oxidizer that fails to size is REPORTED WITH THE REASON, never dropped in
% silence.

clear;
close all;
clc;

%% Parameters and configurations
% This file sits at Optimization/, so the repository root is two levels up.
% Everything below is reached through the path added here.
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(repoRoot));

params = combustion_params();
C = optimizationConstraints();

% ------------ CONFIGURATION -----------

% OX_list = ["O2(L)", "N2O", "N2O4", "IRFNA", "H2O2_90"];
OX_list = ["O2(L)"];

reportFile = fullfile(repoRoot, "docs", "RISULTATI_OTTIMIZZAZIONE.md");

plot_winner = true;       % figures for the winning design
validate_winner = true;   % re-run it through run_combustion_case
                          
grid_search = 350;        % [-] MDF divisions during the search
grid_fine = 700;          % [-] MDF divisions for the ranking

optsA = struct("do_plots", false, "quiet", false, "write_csv", true, ...
    "grid_search", grid_search, "grid_fine", grid_fine);

% Carried in params for the report
params.phaseA_grids = struct("grid_search", grid_search, "grid_fine", grid_fine);
optsB = struct("do_plots", false, "quiet", false, "validate", false);

% Every constraint of the run, printed before anything else
printConstraints(C, params, ...
    struct("grid_search", grid_search, "grid_fine", grid_fine));

%% Enumerate the oxidizers over the whole A -> B chain
n_OX = numel(OX_list);
res = repmat(pack_empty(), 1, n_OX);
run_timer = tic;

for i = 1:n_OX
    fprintf("\n%s\n", repmat('=', 1, 78));
    fprintf(" OXIDIZER %d/%d: %s\n", i, n_OX, OX_list(i));
    fprintf("%s\n", repmat('=', 1, 78));

    try
        thermo = get_thermo(OX_list(i), repoRoot);
    catch ME
        res(i) = pack(OX_list(i), struct([]), ...
            struct("ok", false, "failTag", "thermo: " + string(ME.message)));
        fprintf(" %s: %s\n", OX_list(i), res(i).failTag);
        continue
    end

    shapes = phaseA(thermo, params, C, optsA);
    design = phaseB(shapes, thermo, params, C, optsB);
    res(i) = pack(OX_list(i), shapes, design);
end

fprintf("\nTotal elapsed: %.1f min\n", toc(run_timer)/60);

%% Comparative table
printComparison(res);

%% Best: highest Isp_load among those passing EVERY gate, C13 included
ok = [res.ok];
if ~any(ok)
    error(['No oxidizer produced an admissible design. Read the failTag column ' ...
        'above against the constraint table printed at the top.']);
end
Isp_all = [res.Isp_load];
Isp_all(~ok) = -inf;
[~, ibest] = max(Isp_all);
best = res(ibest);

fprintf("\n%s\n", repmat('=', 1, 78));
fprintf(" BEST: %s, Isp on loaded mass = %.2f s\n", best.oxidizer, best.Isp_load);
fprintf("%s\n", repmat('=', 1, 78));

%% The winner: characteristics, required plots, independent re-run
% reportDesign works on a design that is already sized, so nothing is
% recomputed here beyond the validation run itself.
best.design = reportDesign(best.design, params, ...
    struct("quiet", false, "do_plots", plot_winner, "validate", validate_winner));
res(ibest).design = best.design;

%% Report for the write-up
writeReport(reportFile, res, ibest, C, params, OX_list);
fprintf("\nReport written to %s\n", reportFile);

%% FUNCTIONS

function e = pack_empty()
    % pack_empty
    % Prototype of one result record, so the array can be preallocated.
    % INPUT
    %   None
    % OUTPUT
    %   e / struct / 1x1   Empty result record

    e = struct("oxidizer", "", "ok", false, "failTag", "", "Isp_load", NaN, ...
        "mdot_ox", NaN, "OF_med", NaN, "R_c", NaN, "L", NaN, "LD", NaN, ...
        "Gox0", NaN, "Gox_end", NaN, "sigma", NaN, "drift", NaN, ...
        "geometry", "", "N", NaN, "x1", NaN, "h", NaN, "x2", NaN, ...
        "burn_time", NaN, "thrust0", NaN, "I_tot", NaN, "m_load", NaN, ...
        "shapes", struct([]), "design", struct([]));
end

function e = pack(OX, shapes, design)
    % pack
    % Collect one oxidizer's result into a flat record for the table and the
    % report. A failure keeps its reason.
    % INPUT
    %   OX      / string / 1x1   Oxidizer name
    %   shapes  / struct / 1xK   Phase A ranking (may be empty)
    %   design  / struct / 1x1   Phase B outcome
    % OUTPUT
    %   e       / struct / 1x1   Result record

    e = pack_empty();
    e.oxidizer = string(OX);
    e.shapes = shapes;
    e.design = design;

    if ~isfield(design, "ok") || ~design.ok
        e.ok = false;
        if isfield(design, "failTag")
            e.failTag = string(design.failTag);
        else
            e.failTag = "unknown";
        end
        return
    end

    S = design.S;
    sh = design.shape;
    e.ok = true;
    e.failTag = "OK";
    e.Isp_load = S.Isp_load;
    e.mdot_ox = S.mdot_ox;
    e.OF_med = S.OF_med;
    e.R_c = S.R_c;
    e.L = S.L;
    e.LD = S.L / (2*S.R_c);
    e.Gox0 = S.Gox0;
    e.Gox_end = S.Gox_end;
    e.sigma = S.sigma;
    e.drift = S.drift;
    e.geometry = string(sh.geometry);
    e.N = sh.N;
    e.x1 = sh.x1;
    e.h = sh.h;
    e.x2 = sh.x2;
    e.burn_time = S.burn_time;
    e.thrust0 = S.thrust0;
    e.I_tot = S.I_tot;
    e.m_load = S.m_load;
end

