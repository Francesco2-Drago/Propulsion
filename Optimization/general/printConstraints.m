function printConstraints(C, params, opts)
    % printConstraints
    % Render the whole constraint table of optimizationConstraints at the start
    % of a run. This print plus that one file are meant to be enough for a
    % first-time reader to understand the entire structure of the optimization:
    % what is bounded, what is solved for, what is only watched.
    % INPUT
    %   C      / struct / 1xM   Constraint table. Omit to fetch it here
    %   params / struct / 1x1   Optional. Combustion parameters, used only for
    %                           the values that depend on the discretization
    %                           (C2) or on the fuel (the C4 exponents)
    %   opts   / struct / 1x1   Optional. Fields:
    %                             rationale / logical / print the why-column
    %                                         (default true)
    %                             phase     / string  / "A", "B" or "all"
    %                                         (default "all")
    % OUTPUT
    %   None (prints to stdout)

    if nargin < 1 || isempty(C)
        C = optimizationConstraints();
    end
    if nargin < 2
        params = [];
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, "rationale"), opts.rationale = true; end
    if ~isfield(opts, "phase"), opts.phase = "all"; end
    % Phase A searches on one grid and ranks on another; C2's floor depends on
    % which, so both have to be declared here
    if ~isfield(opts, "grid_search") || isempty(opts.grid_search)
        opts.grid_search = grid_default(params);
    end
    if ~isfield(opts, "grid_fine") || isempty(opts.grid_fine)
        opts.grid_fine = opts.grid_search;
    end

    K = constraintsById(C);
    [~, ~, box] = optimizationConstraints();

    width = 78;
    fprintf("\n%s\n", repmat('=', 1, width));
    fprintf(" OPTIMIZATION CONSTRAINTS - single source: optimizationConstraints.m\n");
    fprintf("%s\n", repmat('=', 1, width));
    fprintf(" No weighted penalties anywhere: every constraint is a bound or a\n");
    fprintf(" boolean gate. Phase A is scale-free (R_c = 1), phase B is the sized\n");
    fprintf(" engine. 'equation' means solved for, not searched.\n");
    fprintf("%s\n", repmat('-', 1, width));
    fprintf(" %-4s %-4s %-9s %-24s %s\n", "id", "ph", "kind", "quantity", "limit");
    fprintf("%s\n", repmat('-', 1, width));

    % Constraint rows, in the binding order C1..C14
    for i = 1:numel(C)
        c = C(i);
        if opts.phase ~= "all" && ~contains(c.phase, opts.phase)
            continue
        end
        for j = 1:numel(c.quantity)
            if j == 1
                fprintf(" %-4s %-4s %-9s %-24s %s\n", c.id, c.phase, c.kind, ...
                    c.quantity(j), limit_string(c, j));
            else
                fprintf(" %-4s %-4s %-9s %-24s %s\n", "", "", "", ...
                    c.quantity(j), limit_string(c, j));
            end
        end
        if opts.rationale
            print_wrapped(c.name + ": " + c.rationale, 6, width);
            if strlength(c.interaction) > 0
                print_wrapped("-> " + c.interaction, 6, width);
            end
        end
    end
    fprintf("%s\n", repmat('-', 1, width));

    % Values that depend on the current discretization or on the fuel
    if ~isempty(params)
        n_rf = params.fuel.n_rf;
        fprintf(" DERIVED AT THE CURRENT SETTINGS\n");

        % C2 depends on the grid, and phase A uses TWO of them: a coarse one to
        % search and a fine one to rank. Both are printed, because the floor
        % halves between them and the optimum sits exactly on whichever is in
        % force. Printing only one hid a factor-two discrepancy between the
        % constraint the search obeyed and the constraint the report declared.
        grids = unique([opts.grid_search, opts.grid_fine], "stable");
        labels = ["search", "ranking"];
        for g = 1:numel(grids)
            if isscalar(grids)
                tag = "search and ranking";
            else
                tag = labels(g);
            end
            fprintf("   C2  grid_divisions = %4d (%-18s) ->  h >= %.5f R_c\n", ...
                grids(g), tag, K.C2.lo / grids(g));
        end
        if ~isscalar(grids)
            fprintf("       the floor MOVES between the two: the shapes are re-optimized\n");
            fprintf("       on the ranking grid, not just re-scored\n");
        end
        fprintf("   C4  G_ox(0) in [%g, %g] %s, i.e. C10 = [%g, %g] with a 5%% margin\n", ...
            K.C4.lo, K.C4.hi, K.C4.units, K.C10.lo, K.C10.hi);
        fprintf("       evaluated on mdot_ox = F/mean[(1+1/OF)*c_eff], exponents %.1f and %.1f\n", ...
            1/(2*n_rf + 1), 2/(2*n_rf + 1));
        fprintf("       at mdot_ox = 10 kg/s this is Gtilde in [%.2f, %.2f]\n", ...
            K.C10.lo / gox_scale(10, n_rf, params.fuel.a_rf, K.C6.lo), ...
            K.C10.hi / gox_scale(10, n_rf, params.fuel.a_rf, K.C6.lo));
        fprintf("   C7  the 50 kN are read as the %s thrust\n", upper(K.C7.value));
        fprintf("%s\n", repmat('-', 1, width));
    end

    % Search box of the decision variables
    fprintf(" SEARCH BOX\n");
    % The floor on h is C2's, not C1's: C1 only says h >= 0, and C2 is always
    % the tighter of the two
    fprintf("   phase A   x1 = ri/R_c in [%.2f, %.2f],  N in %d..%d\n", ...
        box.x1(1), box.x1(2), box.N(1), box.N(2));
    fprintf("             h >= %.5f (C2 on the search grid), C1's h >= 0 never binds\n", ...
        K.C2.lo / opts.grid_search);
    fprintf("   phase B   no box on mdot_ox: it is COMPUTED in phase A from C7,\n");
    fprintf("             and only refined within +/-%.0f%% of that value\n", ...
        100*box.mdot_bracket);
    fprintf("%s\n", repmat('-', 1, width));

    % Diagnostics: computed and reported, deliberately NOT constraints
    fprintf(" DIAGNOSTICS - computed and reported, NOT constraints:\n");
    fprintf("   G_ox(150 s)   unreachable: >= 170 would need a ~45 m grain\n");
    fprintf("   G_ox(end)     follows from mass conservation, almost only L and O/F\n");
    fprintf("   O/F drift     a symptom; its cost is the Isp_med the objective measures\n");
    fprintf("   sliver sigma  enters the objective directly, through the loaded mass\n");
    fprintf("   Phi(0)        the OLD phase A objective: minimizing it maximizes the\n");
    fprintf("                 perimeter, shortens the grain and increases the drift\n");
    fprintf("%s\n\n", repmat('=', 1, width));
end

%% FUNCTIONS

function g = grid_default(params)
    % grid_default
    % Fall back to the grid declared in params when the caller did not say which
    % grids phase A will use.
    % INPUT
    %   params / struct / 1x1   Combustion parameters, possibly empty
    % OUTPUT
    %   g      / double / 1x1   Grid divisions [-]

    if isempty(params) || ~isfield(params, "mdf")
        g = 350;
    else
        g = params.mdf.grid_divisions;
    end
end

function C_scale = gox_scale(mdot_ox, n, a, t_b)
    % gox_scale
    % Scale factor of the initial-flux factorization, Gox0 = C(mdot_ox)*Gtilde.
    % Only used to show what C4 means in shape terms at a representative flow.
    % INPUT
    %   mdot_ox / double / 1x1   Oxidizer mass flow [kg/s]
    %   n       / double / 1x1   Regression exponent [-]
    %   a       / double / 1x1   Regression coefficient [m/s]/[kg/(m2 s)]^n
    %   t_b     / double / 1x1   Burn time [s]
    % OUTPUT
    %   C_scale / double / 1x1   Scale factor [kg/(m2 s)]

    C_scale = mdot_ox^(1/(2*n + 1)) / (t_b*a)^(2/(2*n + 1));
end

function s = limit_string(c, j)
    % limit_string
    % Format one bound pair of a constraint row for the table.
    % INPUT
    %   c / struct / 1x1   Constraint row
    %   j / double / 1x1   Index of the bound pair inside the row
    % OUTPUT
    %   s / string / 1x1   Human-readable limit, units included

    if c.kind == "fixed"
        s = sprintf("= %s", c.value);
        return
    end
    if c.kind == "removed"
        s = "WITHDRAWN";
        return
    end

    lo = c.lo(j);
    hi = c.hi(j);
    u = c.units;
    if u == "-"
        u = "";
    else
        u = " " + u;
    end

    if any(isnan([lo, hi]))
        s = "(derived)";
    elseif c.kind == "equation"
        s = sprintf("= %g%s", lo, u);
    elseif isfinite(lo) && isfinite(hi)
        s = sprintf("in [%g, %g]%s", lo, hi, u);
    elseif isfinite(lo)
        s = sprintf(">= %g%s", lo, u);
    elseif isfinite(hi)
        s = sprintf("<= %g%s", hi, u);
    else
        s = "unbounded";
    end
end

function print_wrapped(text, indent, width)
    % print_wrapped
    % Print a rationale wrapped to the table width, so the why-column stays
    % readable next to the numbers.
    % INPUT
    %   text   / string / 1x1   Text to wrap
    %   indent / double / 1x1   Left indentation [characters]
    %   width  / double / 1x1   Total line width [characters]
    % OUTPUT
    %   None (prints to stdout)

    words = split(strtrim(text));
    pad = repmat(' ', 1, indent);
    line = "";
    for i = 1:numel(words)
        if strlength(line) == 0
            candidate = words(i);
        else
            candidate = line + " " + words(i);
        end
        if strlength(candidate) + indent > width
            fprintf("%s%s\n", pad, line);
            line = words(i);
        else
            line = candidate;
        end
    end
    if strlength(line) > 0
        fprintf("%s%s\n", pad, line);
    end
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
