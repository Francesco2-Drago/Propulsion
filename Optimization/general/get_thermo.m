function thermo = get_thermo(OX, repoRoot)
    % get_thermo
    % Load the CEA interpolants of one oxidizer and wrap them as handles of
    % (O/F, p_c[bar]). The oxidizer parametrizes the whole A -> B chain, so the
    % thermochemistry is loaded once per oxidizer and passed down into both
    % phases rather than being reloaded inside them.
    %
    % The lookup DOMAIN is not attached here: it comes from the constraint
    % table (C5) and from the pressure extent of the tables, and is added by
    % engineContext, which owns both.
    % INPUT
    %   OX       / string / 1x1   Oxidizer name, matching the LookupTable folder
    %   repoRoot / string / 1x1   Repository root
    % OUTPUT
    %   thermo   / struct / 1x1   name, Tc_fun, R_fun, k_fun [handles of
    %                             (of, p_bar)] and the R*Tc derivatives the
    %                             fine ODE needs

    OX = string(OX);
    f = fullfile(repoRoot, "LookupTable", OX, "interpolant_" + OX + ".mat");
    if ~isfile(f)
        error("get_thermo:missingLookup", ...
            "No interpolant for oxidizer '%s' at %s.", OX, f);
    end
    d = load(f);

    T_interp = d.interpolant.T_interp;             % [K]
    R_interp = d.interpolant.R_interp;             % [J/kgK]
    g_interp = d.interpolant.gamma_interp;         % [-]

    thermo.name = OX;
    thermo.Tc_fun = @(of, pb) T_interp(pb, of);
    thermo.R_fun = @(of, pb) R_interp(pb, of);
    thermo.k_fun = @(of, pb) g_interp(pb, of);

    % Derivatives of R*Tc, for the fine ODE of the chamber-pressure model
    dp = 1e-3; dof = 1e-3;
    RT = @(of, pb) R_interp(pb, of) * T_interp(pb, of);
    thermo.dRTdp_fun = @(of, pb) (RT(of, pb + dp) - RT(of, pb - dp)) / (2*dp);
    thermo.dRTdOF_fun = @(of, pb) (RT(of + dof, pb) - RT(of - dof, pb)) / (2*dof);
end
