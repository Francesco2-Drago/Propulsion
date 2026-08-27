function design = reportDesign(design, params, opts)
    % reportDesign
    % Present one sized engine: the characteristics table, the diagnostic
    % figures, and the independent re-run through run_combustion_case that
    % produces the histories the assignment asks for.
    %
    % Kept separate from phaseB so it can be called on a design that has
    % ALREADY been sized, without paying for the sizing again. main_optimization
    % uses that to plot the winner only, instead of every oxidizer it tried.
    %
    % The three options are INDEPENDENT: printing, plotting and validating are
    % three different things and none of them gates the others.
    % INPUT
    %   design / struct / 1x1   Output of phaseB: S, shape, ctx, sweep,
    %                           mdot_grid
    %   params / struct / 1x1   Combustion parameters
    %   opts   / struct / 1x1   Optional: quiet (default false), do_plots
    %                           (default true), validate (default false)
    % OUTPUT
    %   design / struct / 1x1   The same design, with .validation added when
    %                           the validation run was requested

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, "quiet") || isempty(opts.quiet), opts.quiet = false; end
    if ~isfield(opts, "do_plots") || isempty(opts.do_plots), opts.do_plots = true; end
    if ~isfield(opts, "validate") || isempty(opts.validate), opts.validate = false; end

    if ~isfield(design, "ok") || ~design.ok
        if ~opts.quiet
            fprintf("No design to report: %s\n", string(design.failTag));
        end
        return
    end

    S = design.S;
    shape = design.shape;
    ctx = design.ctx;

    % ------------ CHARACTERISTICS ------------
    if ~opts.quiet
        print_characteristics(S, shape, ctx);
    end

    % ------------ DIAGNOSTIC FIGURES ------------
    if opts.do_plots
        plot_sweep_landscape(design, S);
        plot_grain(shape, S, ctx);
    end

    % ------------ INDEPENDENT RE-RUN ------------
    if opts.validate
        design = validate_design(design, params, S, shape, ctx, opts);
    end
end

%% FUNCTIONS

function print_characteristics(S, shape, ctx)
    % print_characteristics
    % The full engine table.
    % INPUT
    %   S     / struct / 1x1   Sized engine
    %   shape / struct / 1x1   Shape record
    %   ctx   / struct / 1x1   Evaluation context
    % OUTPUT
    %   None (prints to stdout)

    fprintf("\n==================== OPTIMAL ENGINE ====================\n");
    fprintf(" Oxidizer                 : %s\n", ctx.oxidizer);
    fprintf(" Grain geometry           : %s\n", shape.geometry);
    fprintf(" Objective, Isp on loaded : %.2f s\n", S.Isp_load);
    fprintf(" Isp on burnt mass        : %.2f s\n", S.Isp_med);
    fprintf("\n --- Operating point ---\n");
    fprintf(" mdot_ox (const)          : %.3f kg/s\n", S.mdot_ox);
    fprintf(" O/F  (start -> end)      : %.2f -> %.2f\n", S.of0, S.of_end);
    fprintf(" O/F  mean (mass)         : %.3f\n", S.OF_med);
    fprintf(" G_ox (start -> end)      : %.0f -> %.0f kg/m2 s\n", S.Gox0, S.Gox_end);
    fprintf(" peak chamber pressure    : %.2f bar\n", S.p_peak*1e-5);
    fprintf("\n --- Geometry ---\n");
    fprintf(" casing diameter (2 R_c)  : %.3f m\n", 2*S.R_c);
    fprintf(" grain length L           : %.3f m\n", S.L);
    fprintf(" L/D                      : %.2f\n", S.L/(2*S.R_c));
    fprintf(" port inner diameter      : %.3f m\n", 2*shape.x1*S.R_c);
    if shape.geometry ~= "cylinder"
        fprintf(" port tip diameter        : %.3f m\n", 2*shape.x2*S.R_c);
        fprintf(" number of tips           : %d\n", shape.N);
    end
    fprintf(" throat diameter          : %.4f m  (At = %.4e m2)\n", ...
        S.throat_diameter, S.At);
    fprintf(" nozzle area ratio        : %d\n", ctx.engine.eps);
    fprintf(" residual web (C12)       : %.1f mm\n", S.web_residual*1e3);
    fprintf("\n --- Masses ---\n");
    fprintf(" oxidizer loaded          : %.1f kg\n", S.mdot_ox*S.burn_time);
    fprintf(" fuel loaded              : %.1f kg\n", S.m_load - S.mdot_ox*S.burn_time);
    fprintf(" LOADED mass              : %.1f kg\n", S.m_load);
    fprintf(" mass through the nozzle  : %.1f kg\n", S.m_prop);
    fprintf(" sliver (unburnt fuel)    : %.2f %%\n", 100*S.sigma);
    fprintf("\n --- Performance ---\n");
    fprintf(" initial thrust           : %.2f kN\n", S.thrust0*1e-3);
    fprintf(" mean thrust              : %.2f kN\n", S.mean_thrust*1e-3);
    fprintf(" total impulse            : %.3f MNs\n", S.I_tot*1e-6);
    fprintf(" burn time                : %.1f s\n", S.burn_time);
    fprintf("\n --- Diagnostics (NOT constraints) ---\n");
    fprintf(" G_ox at 150 s            : %.1f kg/m2 s\n", S.Gox_150);
    fprintf(" G_ox minimum             : %.1f kg/m2 s\n", S.Gox_min);
    fprintf(" O/F drift (max/min)      : %.2f\n", S.drift);
    fprintf(" regression rate rf       : %.3f -> %.3f mm/s\n", S.rf0*1e3, S.rf_end*1e3);
    fprintf("========================================================\n");
end

function plot_sweep_landscape(design, S)
    % plot_sweep_landscape
    % The objective against the oxidizer flow, with the points the gates
    % refused shown in grey: the shape of the admissible window is worth seeing.
    % INPUT
    %   design / struct / 1x1   Design record with sweep and mdot_grid
    %   S      / struct / 1x1   Sized engine
    % OUTPUT
    %   None (creates a figure)

    if ~isfield(design, "sweep") || isempty(design.sweep)
        return
    end

    figure("Name", "Phase B - Isp on loaded mass vs mdot_ox", "Color", "w");
    hold on
    Isp_plot = [design.sweep.Isp_load];
    feas_plot = [design.sweep.feasible];
    if any(~feas_plot)
        plot(design.mdot_grid(~feas_plot), Isp_plot(~feas_plot), "o", ...
            "Color", [0.6 0.6 0.6], "DisplayName", "gate failed");
    end
    if any(feas_plot)
        plot(design.mdot_grid(feas_plot), Isp_plot(feas_plot), "o-", ...
            "LineWidth", 1.5, "Color", [0 0.45 0.74], "DisplayName", "admissible");
    end
    plot(S.mdot_ox, S.Isp_load, "p", "MarkerSize", 15, ...
        "MarkerFaceColor", [0.85 0.33 0.10], "MarkerEdgeColor", "k", ...
        "DisplayName", "optimum");
    grid on
    xlabel("$\dot m_{ox}$ [kg/s]", "Interpreter", "latex");
    ylabel("$I_{sp}$ on loaded mass [s]", "Interpreter", "latex");
    title("Phase B objective vs oxidizer flow");
    legend("Location", "best");
end

function plot_grain(shape, S, ctx)
    % plot_grain
    % The grain cross-section at the sized casing radius.
    % INPUT
    %   shape / struct / 1x1   Shape record
    %   S     / struct / 1x1   Sized engine
    %   ctx   / struct / 1x1   Evaluation context
    % OUTPUT
    %   None (creates a figure)

    figure("Name", "Phase B - optimal engine grain", "Color", "w");
    th = linspace(0, 2*pi, 400);
    if shape.geometry == "cylinder"
        r_port = shape.x1 * S.R_c;
        port = [r_port*cos(th).', r_port*sin(th).'];
    else
        md.inner_diameter = 2 * shape.x1 * S.R_c;
        md.outer_diameter = 2 * shape.x2 * S.R_c;
        md.n_tips = shape.N;
        half = make_mesh0("star", md, 800, "cartesian");
        port = [half; flipud([half(:,1), -half(:,2)])];
    end
    fill(S.R_c*cos(th), S.R_c*sin(th), [0.82 0.72 0.47], ...
        "EdgeColor", "k", "LineWidth", 1.3);
    hold on
    fill(port(:,1), port(:,2), "w", "EdgeColor", [0.20 0.20 0.20], "LineWidth", 1.1);
    axis equal
    axis(S.R_c * [-1.1 1.1 -1.1 1.1])
    set(gca, "XTick", [], "YTick", [])
    title(sprintf("%s | 2R_c = %.2f m, L = %.2f m", ctx.oxidizer, 2*S.R_c, S.L), ...
        "Interpreter", "none");
end

function design = validate_design(design, params, S, shape, ctx, opts)
    % validate_design
    % Re-run the winning design through run_combustion_case, which rebuilds the
    % MDF from scratch, and compare. This is the coherence check of the brief
    % 7.3, and it is also what produces the histories the assignment requires:
    % chamber pressure and temperature, O/F and G_ox, thrust and specific
    % impulse, plus the regression rate.
    % INPUT
    %   design / struct / 1x1   Design record
    %   params / struct / 1x1   Combustion parameters
    %   S      / struct / 1x1   Sized engine
    %   shape  / struct / 1x1   Shape record
    %   ctx    / struct / 1x1   Evaluation context
    %   opts   / struct / 1x1   Options
    % OUTPUT
    %   design / struct / 1x1   With .validation added

    params_win = params;
    params_win.combustion.oxidizer_lookup = ctx.oxidizer;
    params_win.combustion.mdot_ox = S.mdot_ox;
    params_win.engine.ext_diameter = 2 * S.R_c;
    params_win.engine.throat_diameter = S.throat_diameter;
    params_win.engine.chamber_length = S.L;
    if shape.geometry == "cylinder"
        params_win.geometry.type = "cylinder";
        params_win.geometry.meshdata = struct("diameter", 2 * shape.x1 * S.R_c);
    else
        params_win.geometry.type = "star";
        params_win.geometry.meshdata = struct( ...
            "inner_diameter", 2 * shape.x1 * S.R_c, ...
            "outer_diameter", 2 * shape.x2 * S.R_c, ...
            "n_tips", shape.N);
    end
    % The rebuild has to stop where the design stops. run_combustion_case takes
    % b_max straight from the MDF and knows nothing about C12, so left alone it
    % burns all the way to the casing, where the perimeter collapses and the
    % O/F leaves the CEA table: it would crash in the validating call of its
    % STEP 4, and even if it did not it would be integrating a different
    % problem. Capping the horizon at the design's own burn time makes the two
    % cover the same interval. The comparison of the burn time below is still
    % informative: if the rebuild hits an event of its own it stops EARLIER
    % than the cap, and the delta shows it.
    params_win.time.tmax = S.burn_time;
    params_win.time.dt_output = 0.25;

    if ~opts.quiet
        fprintf("\nRe-running the design through run_combustion_case (MDF rebuilt from scratch)...\n");
    end
    [t_val, in_time_val, summary_val] = run_combustion_case(params_win, ...
        struct("do_plots", opts.do_plots, "do_animation", false, ...
        "quiet", opts.quiet));

    if ~opts.quiet
        fprintf("\n --- Sizing model vs full rebuild (should agree within 1%%) ---\n");
        fprintf(" %-18s %12s %12s %8s\n", "quantity", "sizing", "rebuild", "delta");
        compare_row("burn time [s]", S.burn_time, summary_val.burn_time);
        compare_row("mean thrust [kN]", S.mean_thrust*1e-3, summary_val.mean_thrust*1e-3);
        compare_row("peak p_c [bar]", S.p_peak*1e-5, summary_val.max_pressure*1e-5);
        compare_row("total impulse [MNs]", S.I_tot*1e-6, summary_val.total_impulse*1e-6);
    end

    design.validation.time = t_val;
    design.validation.in_time = in_time_val;
    design.validation.summary = summary_val;

    % Regression rate history, explicitly requested by the assignment and not
    % among the plots run_combustion_case already draws
    if opts.do_plots
        figure("Name", "Regression rate", "Color", "w");
        plot(t_val, params.fuel.a_rf*in_time_val.GOX.^params.fuel.n_rf*1e3, ...
            "LineWidth", 1.6, "Color", [0.85 0.33 0.10]);
        grid on
        xlabel("time, t [s]");
        ylabel("r_f [mm/s]");
        title("Fuel regression rate");
    end
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
