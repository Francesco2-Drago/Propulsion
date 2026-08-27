function lookup = cylinderLookup(x1, params)
    % cylinderLookup
    % Closed-form burnback lookup of a circular port, in the same format
    % build_shape_lookup returns, so the cylinder can be scored by the very
    % same shapeMerit as the star family and compete with it on equal terms.
    %
    % A concentric circular port is the one shape with no pathological tail:
    % every point of the surface reaches the casing at the same instant, so
    % Ap and Pb grow together, the perimeter ends at 2*pi*R_c instead of
    % collapsing to zero, and the O/F never runs off the CEA table. The star
    % family cannot reach this configuration, because C2 needs a resolvable
    % tip; the optimizer pressing against that floor is the symptom.
    %
    % Under uniform normal regression:
    %   r(b) = r0 + b,   Ap = pi*r^2,   Pb = 2*pi*r,   dPb/db = 2*pi
    % the last identity being Steiner's formula for a convex port, which the
    % acceptance tests check against the MDF.
    % INPUT
    %   x1     / double / 1x1   Port radius ratio r0/R_c [-]
    %   params / struct / 1x1   Combustion parameters: engine.ext_diameter [m]
    %                           sets the casing radius, mdf.n_lookup the number
    %                           of burnback levels
    % OUTPUT
    %   lookup / struct / 1x1   b [m], Ap [m2], perim [m], Ap_raw, dAp_db,
    %                           dperim_db, b_max [m]

    R_c = params.engine.ext_diameter / 2;      % [m]
    if x1 <= 0 || x1 >= 1
        error("cylinderLookup:badPort", ...
            "The port ratio must lie strictly between 0 and 1, got %g.", x1);
    end

    if isfield(params, "mdf") && isfield(params.mdf, "n_lookup")
        n_pts = params.mdf.n_lookup;
    else
        n_pts = 601;
    end

    r0 = x1 * R_c;                             % [m] initial port radius
    b = linspace(0, R_c - r0, n_pts).';        % [m] burnback, up to the casing
    r = r0 + b;                                % [m] instantaneous port radius

    lookup.b = b;
    lookup.Ap = pi * r.^2;                     % [m2]
    lookup.perim = 2*pi * r;                   % [m]
    lookup.Ap_raw = lookup.Ap;
    lookup.dAp_db = lookup.perim;              % dAp/db = perimeter
    lookup.dperim_db = 2*pi * ones(size(b));   % [-]  Steiner
    lookup.b_max = b(end);                     % [m]
end
