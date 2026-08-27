% PLOT_RESULTS - figures, validation and report from results ALREADY in the
% workspace, without re-running the optimization.
%
% Run this after main_optimization (or after a bare phaseA/phaseB pair) when
% you want the plots, the independent re-run, or a fresh report, and you do not
% want to pay for the search again. Everything here works on a design that has
% already been sized: the only thing that costs time is the validation re-run,
% about a minute, and it is optional.
%
% It looks for, in order:
%   best.design       what main_optimization leaves behind
%   res(ibest).design idem, if best was cleared
%   design            what a bare phaseB call leaves behind
%
% Usage:
%   plotResults                     % figures + validation + report
%   do_validate = false; plotResults   % skip the re-run, figures only

%% What is available in the workspace
if ~exist("params", "var")
    fprintf("params not in the workspace, rebuilding from combustion_params().\n");
    params = combustion_params();
end
if ~exist("C", "var")
    C = optimizationConstraints();
end

% ------------ CONFIGURATION -----------
% Set these before calling if you want something other than the default
if ~exist("do_plots", "var"),    do_plots = true;    end
if ~exist("do_validate", "var"), do_validate = true; end
if ~exist("do_report", "var"),   do_report = true;   end

%% Locate the design
design_found = [];
design_source = "";

if exist("best", "var") && isstruct(best) && isfield(best, "design") ...
        && isfield(best.design, "ok") && best.design.ok
    design_found = best.design;
    design_source = "best.design";
elseif exist("res", "var") && isstruct(res) && exist("ibest", "var") ...
        && isfield(res(ibest), "design") && isfield(res(ibest).design, "ok") ...
        && res(ibest).design.ok
    design_found = res(ibest).design;
    design_source = sprintf("res(%d).design", ibest);
elseif exist("design", "var") && isstruct(design) && isfield(design, "ok") ...
        && design.ok
    design_found = design;
    design_source = "design";
end

if isempty(design_found)
    error(['No sized design found in the workspace. Expected one of: best.design, ' ...
        'res(ibest).design, or design. Run main_optimization first, or phaseA ' ...
        'followed by phaseB.']);
end

fprintf("\nUsing %s: %s, %s, Isp_load = %.2f s\n", design_source, ...
    design_found.ctx.oxidizer, design_found.shape.geometry, ...
    design_found.S.Isp_load);

%% Comparative table, if the full enumeration is available
if exist("res", "var") && isstruct(res) && isfield(res, "oxidizer")
    printComparison(res);
end

%% Characteristics, figures, and the independent re-run
% The three switches are independent: none of them gates the others.
design_found = reportDesign(design_found, params, struct( ...
    "quiet", false, "do_plots", do_plots, "validate", do_validate));

% Put the enriched design back where it came from, so the validation histories
% stay available in the workspace
if design_source == "best.design"
    best.design = design_found;
    if exist("res", "var") && exist("ibest", "var")
        res(ibest).design = design_found;
    end
elseif startsWith(design_source, "res(")
    res(ibest).design = design_found;
else
    design = design_found;
end

%% Report for the write-up
if do_report
    if ~exist("res", "var") || ~isstruct(res) || ~isfield(res, "oxidizer")
        fprintf(["\nSkipping the report: it needs the full res array from " ...
            "main_optimization, not just one design.\n"]);
    else
        if ~exist("ibest", "var")
            Isp_all = [res.Isp_load];
            Isp_all(~[res.ok]) = -inf;
            [~, ibest] = max(Isp_all);
        end
        if ~exist("OX_list", "var")
            OX_list = [res.oxidizer];
        end
        if ~exist("reportFile", "var")
            reportFile = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), ...
                "docs", "RISULTATI_OTTIMIZZAZIONE.md");
        end
        writeReport(reportFile, res, ibest, C, params, OX_list);
        fprintf("\nReport written to %s\n", reportFile);
    end
end

fprintf("\nDone. %d figures open.\n", numel(findobj("Type", "figure")));
