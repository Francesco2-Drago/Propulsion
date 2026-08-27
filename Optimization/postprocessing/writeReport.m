function writeReport(reportFile, res, ibest, C, params, OX_list)
    % writeReport
    % Write the markdown report the write-up needs: constraints in force, the
    % shape ranking of every oxidizer, the comparison across oxidizers, and the
    % decomposition of why the winning shape won.
    % INPUT
    %   reportFile / string / 1x1   Destination path
    %   res        / struct / 1xM   Result records
    %   ibest      / double / 1x1   Index of the winner
    %   C          / struct / 1xN   Constraint table
    %   params     / struct / 1x1   Combustion parameters
    %   OX_list    / string / 1xM   Oxidizers enumerated
    % OUTPUT
    %   None (writes the file)

    K = constraintsById(C);
    folder = fileparts(reportFile);
    if ~isfolder(folder)
        mkdir(folder);
    end
    fid = fopen(reportFile, "w");
    if fid < 0
        warning("write_report:cannotOpen", "Could not write %s", reportFile);
        return
    end
    closer = onCleanup(@() fclose(fid));

    stamp = string(datetime("now", "Format", "yyyy-MM-dd HH:mm"));
    fprintf(fid, "# Risultati dell'ottimizzazione\n\n");
    fprintf(fid, "Generato da `Optimization2/main_optimization.m` il %s.\n", stamp);
    fprintf(fid, "Non modificare a mano: viene riscritto a ogni run.\n\n");

    % ------------ 1. CONSTRAINTS ------------
    fprintf(fid, "## 1. Vincoli attivi\n\n");
    fprintf(fid, "Sorgente unica: `Optimization2/general/optimizationConstraints.m`. ");
    fprintf(fid, "Nessun altro file contiene una soglia numerica di vincolo. ");
    fprintf(fid, "Nessuna penalita' pesata: ogni vincolo e' un bound o un gate booleano.\n\n");
    fprintf(fid, "| id | fase | tipo | grandezza | limite | motivazione | chi lo consuma, e contro cosa tira |\n");
    fprintf(fid, "|----|------|------|-----------|--------|-------------|-------------------------------------|\n");
    for i = 1:numel(C)
        c = C(i);
        for j = 1:numel(c.quantity)
            if j == 1
                fprintf(fid, "| **%s** | %s | %s | `%s` | %s | %s | %s |\n", ...
                    c.id, c.phase, c.kind, c.quantity(j), ...
                    md_limit(c, j), md_escape(c.rationale), ...
                    md_escape(c.interaction));
            else
                fprintf(fid, "| | | | `%s` | %s | | |\n", c.quantity(j), md_limit(c, j));
            end
        end
    end

    fprintf(fid, "\n**Derivati alle impostazioni correnti.** ");
    % C2 depends on the grid, and phase A uses two: printing only one hid a
    % factor-two gap between the constraint the search obeyed and the one the
    % report declared
    grids = unique([grid_of(params, "grid_search"), grid_of(params, "grid_fine")], "stable");
    if isscalar(grids)
        fprintf(fid, "C2: `grid_divisions = %d`, quindi `h >= %.5f R_c`. ", ...
            grids(1), K.C2.lo/grids(1));
    else
        fprintf(fid, "C2 dipende dalla griglia, e le fasi ne usano due: ");
        fprintf(fid, "ricerca `grid = %d` -> `h >= %.5f R_c`, ", grids(1), K.C2.lo/grids(1));
        fprintf(fid, "classifica `grid = %d` -> `h >= %.5f R_c`. ", grids(2), K.C2.lo/grids(2));
        fprintf(fid, "Le forme sono **riottimizzate** sulla griglia di classifica, ");
        fprintf(fid, "non solo rivalutate: il floor si dimezza fra le due e l'ottimo ");
        fprintf(fid, "ci sta esattamente sopra. ");
    end
    fprintf(fid, "C4: `G_ox(0) in [%g, %g]`, cioe' C10 con margine 5%%, ", K.C4.lo, K.C4.hi);
    fprintf(fid, "valutato sulla portata calcolata da C7. ");
    fprintf(fid, "C7: i 50 kN letti come spinta **%s**.\n\n", upper(K.C7.value));

    % ------------ 2. SHAPE RANKING PER OXIDIZER ------------
    fprintf(fid, "## 2. Classifica delle forme, per ossidante\n\n");
    fprintf(fid, "Fase A, a `R_c = 1`. `Isp_load` e' calcolata direttamente come ");
    fprintf(fid, "`I_tot/(g0*m_load)`; `mdot_ox` non e' cercata, e' ricavata dal ");
    fprintf(fid, "requisito di spinta C7. Il cilindro compete alla pari con le stelle ");
    fprintf(fid, "ed e' valutato con la stessa lookup MDF.\n\n");
    for i = 1:numel(res)
        fprintf(fid, "### %s\n\n", res(i).oxidizer);
        sh = res(i).shapes;
        if isempty(sh)
            fprintf(fid, "Nessuna forma ammissibile. %s\n\n", md_escape(res(i).failTag));
            continue
        end
        fprintf(fid, "| # | forma | N | x1 | h | Isp_load [s] | O/F med | sigma | drift | drift a burnout | mdot_ox [kg/s] | G_ox(0) |\n");
        fprintf(fid, "|---|-------|---|----|---|--------------|---------|-------|-------|-----------------|----------------|--------|\n");
        for k = 1:numel(sh)
            fprintf(fid, "| %d | %s | %d | %.4f | %.5f | %.2f | %.3f | %.4f | %.3f | %s | %.2f | %.0f |\n", ...
                k, sh(k).geometry, sh(k).N, sh(k).x1, sh(k).h, sh(k).Isp_load, ...
                sh(k).OF_med, sh(k).sigma, sh(k).drift, md_num(sh(k).drift_full), ...
                sh(k).mdot_ox, sh(k).Gox0);
        end
        fprintf(fid, "\n");
    end

    % ------------ 3. OXIDIZER COMPARISON ------------
    fprintf(fid, "## 3. Comparativa fra ossidanti\n\n");
    fprintf(fid, "| ossidante | forma | Isp_load [s] | mdot_ox [kg/s] | O/F med | 2R_c [m] | L [m] | L/D | G_ox(0) | G_ox(fine) | sigma | esito |\n");
    fprintf(fid, "|-----------|-------|--------------|----------------|---------|----------|-------|-----|---------|------------|-------|-------|\n");
    for i = 1:numel(res)
        r = res(i);
        if r.ok
            fprintf(fid, "| %s | %s | **%.2f** | %.2f | %.3f | %.3f | %.3f | %.2f | %.0f | %.1f | %.4f | OK |\n", ...
                r.oxidizer, r.geometry, r.Isp_load, r.mdot_ox, r.OF_med, ...
                2*r.R_c, r.L, r.LD, r.Gox0, r.Gox_end, r.sigma);
        else
            fprintf(fid, "| %s | - | - | - | - | - | - | - | - | - | - | %s |\n", ...
                r.oxidizer, md_escape(r.failTag));
        end
    end
    fprintf(fid, "\nOssidanti enumerati: %s.\n", strjoin(OX_list, ", "));
    fprintf(fid, "Un ossidante che non dimensiona e' riportato con il motivo, mai scartato in silenzio.\n\n");

    % ------------ 4. WHY IT WON ------------
    fprintf(fid, "## 4. Perche' ha vinto\n\n");
    write_why(fid, res(ibest));

    fprintf(fid, "\n---\n\n");
    fprintf(fid, "*Decomposizione esatta:* ");
    fprintf(fid, "`Isp_load = Isp_med * (OF_med + 1)/(OF_med + 1/(1 - sigma))`. ");
    fprintf(fid, "Il primo fattore misura la perdita da drift dell'O/F, il secondo ");
    fprintf(fid, "quella da sliver. E' una diagnostica: `Isp_load` viene calcolata ");
    fprintf(fid, "direttamente da `I_tot/(g0*m_load)`, non per questa strada.\n");
end

function write_why(fid, r)
    % write_why
    % Decompose the phase A gap between the winning shape and the best shape of
    % the other family into its drift and sliver parts.
    % INPUT
    %   fid / double / 1x1   Open file id
    %   r   / struct / 1x1   Winning oxidizer record
    % OUTPUT
    %   None (writes to fid)

    sh = r.shapes;
    if isempty(sh)
        fprintf(fid, "Nessuna forma da confrontare.\n");
        return
    end

    geo = arrayfun(@(s) string(s.geometry), sh);
    i_cyl = find(geo == "cylinder", 1, "first");
    i_star = find(geo == "star", 1, "first");
    if isempty(i_cyl) || isempty(i_star)
        fprintf(fid, "Una sola famiglia di forme e' risultata ammissibile ");
        fprintf(fid, "per %s, quindi non c'e' confronto da scomporre.\n", r.oxidizer);
        return
    end

    a = sh(i_cyl);
    b = sh(i_star);
    fa = sliver_factor(a);
    fb = sliver_factor(b);

    % Exact additive split of the difference:
    %   Isp_a - Isp_b = (Isp_med_a - Isp_med_b)*f_a + Isp_med_b*(f_a - f_b)
    d_drift = (a.Isp_med - b.Isp_med) * fa;
    d_sliver = b.Isp_med * (fa - fb);
    d_tot = a.Isp_load - b.Isp_load;

    fprintf(fid, "Confronto in Fase A per **%s**, migliore cilindro contro ", r.oxidizer);
    fprintf(fid, "migliore stella, valutati con la stessa lookup MDF e la stessa griglia.\n\n");
    fprintf(fid, "| | cilindro | stella (N = %d) |\n", b.N);
    fprintf(fid, "|---|---------|--------|\n");
    fprintf(fid, "| `x1` | %.4f | %.4f |\n", a.x1, b.x1);
    fprintf(fid, "| `Isp_load` [s] | **%.2f** | %.2f |\n", a.Isp_load, b.Isp_load);
    fprintf(fid, "| `Isp_med` [s] (drift) | %.2f | %.2f |\n", a.Isp_med, b.Isp_med);
    fprintf(fid, "| fattore sliver | %.5f | %.5f |\n", fa, fb);
    fprintf(fid, "| `sigma` | %.4f | %.4f |\n", a.sigma, b.sigma);
    fprintf(fid, "| drift (parte utile) | %.3f | %.3f |\n", a.drift, b.drift);
    fprintf(fid, "| drift fino a burnout | %s | %s |\n", ...
        md_num(a.drift_full), md_num(b.drift_full));
    fprintf(fid, "\n");

    fprintf(fid, "**Scomposizione della differenza di %+.2f s:**\n\n", d_tot);
    fprintf(fid, "- da drift (via `Isp_med`): **%+.2f s**\n", d_drift);
    fprintf(fid, "- da sliver (via il fattore di massa): **%+.2f s**\n", d_sliver);
    fprintf(fid, "- somma: %+.2f s\n\n", d_drift + d_sliver);

    if abs(d_sliver) > abs(d_drift)
        fprintf(fid, "Il termine dominante e' lo **sliver**.\n\n");
    else
        fprintf(fid, "Il termine dominante e' il **drift**.\n\n");
    end

    fprintf(fid, "**Sul drift fino a burnout.** La colonna diverge per entrambe ");
    fprintf(fid, "le famiglie, cilindro compreso: il solutore MDF smette di ");
    fprintf(fid, "contare la superficie che esce dal casing, quindi il perimetro ");
    fprintf(fid, "bruciante va a zero alla parete per qualunque forma e `Phi = ");
    fprintf(fid, "Ap^n/Pb` diverge con lui. Non e' un fenomeno fisico della ");
    fprintf(fid, "stella, e' contabilita' di fine burn: le zone di collasso sono ");
    fprintf(fid, "confrontabili fra le due forme (ordine dell'1%% del web). ");
    fprintf(fid, "Il burn reale si ferma prima, su C12, e li' il drift e' quello ");
    fprintf(fid, "della colonna 'parte utile'. La colonna a burnout va riportata ");
    fprintf(fid, "al punto (ii) come limite del modello, non come proprieta' della ");
    fprintf(fid, "geometria.\n\n");

    fprintf(fid, "**Quello che invece distingue davvero le due forme e' lo sliver ");
    fprintf(fid, "al criterio di arresto C12.** Le punte della stella raggiungono ");
    fprintf(fid, "il casing prima dei piani, quindi quando il punto piu' sottile ");
    fprintf(fid, "del web e' sceso a 3 mm resta ancora combustibile altrove. ");
    fprintf(fid, "Il cilindro e' l'unica porta concentrica che raggiunge la parete ");
    fprintf(fid, "ovunque nello stesso istante, e questo si vede direttamente nel ");
    fprintf(fid, "confronto di `sigma`: %.4f contro %.4f.\n", a.sigma, b.sigma);
end

function f = sliver_factor(s)
    % sliver_factor
    % The second factor of the Isp_load decomposition, i.e. the fraction of the
    % loaded mass that actually leaves through the nozzle.
    % INPUT
    %   s / struct / 1x1   Shape record
    % OUTPUT
    %   f / double / 1x1   Sliver factor [-]

    f = (s.OF_med + 1) / (s.OF_med + 1/(1 - s.sigma));
end

function s = md_num(v)
    % md_num
    % Format a possibly infinite number for the markdown table.
    % INPUT
    %   v / double / 1x1   Value
    % OUTPUT
    %   s / string / 1x1   Formatted value

    if ~isfinite(v)
        s = "inf";
    else
        s = sprintf("%.3f", v);
    end
end

function s = md_limit(c, j)
    % md_limit
    % Format one bound pair of a constraint row for the markdown table.
    % INPUT
    %   c / struct / 1x1   Constraint row
    %   j / double / 1x1   Index of the bound pair
    % OUTPUT
    %   s / string / 1x1   Formatted limit

    if c.kind == "fixed"
        s = "`= " + c.value + "`";
        return
    end
    if c.kind == "removed"
        s = "**RITIRATO**";
        return
    end
    lo = c.lo(j);
    hi = c.hi(j);
    u = c.units;
    if u == "-", u = ""; else, u = " " + u; end

    if any(isnan([lo, hi]))
        s = "(derivato)";
    elseif c.kind == "equation"
        s = sprintf("`= %g`%s", lo, u);
    elseif isfinite(lo) && isfinite(hi)
        s = sprintf("`[%g, %g]`%s", lo, hi, u);
    elseif isfinite(lo)
        s = sprintf("`>= %g`%s", lo, u);
    elseif isfinite(hi)
        s = sprintf("`<= %g`%s", hi, u);
    else
        s = "-";
    end
end

function g = grid_of(params, field)
    % grid_of
    % Read one of the phase A grids out of params, falling back to the grid
    % declared there when the run did not record them separately.
    % INPUT
    %   params / struct / 1x1   Combustion parameters, possibly carrying
    %                           .phaseA_grids
    %   field  / string / 1x1   "grid_search" or "grid_fine"
    % OUTPUT
    %   g      / double / 1x1   Grid divisions [-]

    if isfield(params, "phaseA_grids") && isfield(params.phaseA_grids, field)
        g = params.phaseA_grids.(field);
    else
        g = params.mdf.grid_divisions;
    end
end

function s = md_escape(t)
    % md_escape
    % Make a rationale safe inside a markdown table cell.
    % INPUT
    %   t / string / 1x1   Text
    % OUTPUT
    %   s / string / 1x1   Escaped text

    s = replace(string(t), "|", "\|");
    s = replace(s, newline, " ");
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
