function lookup = build_shape_lookup(meshdata, params)
    % build_shape_lookup
    % Build the MDF burnback geometry lookup (burnback, port area, burning
    % perimeter) for a given grain cross-section. This is the reusable core
    % shared by the shape optimization and the later length optimization:
    % it does NOT depend on chamber length nor on the oxidizer.
    % INPUT
    %   meshdata / struct / 1x1   Grain cross-section data (inner/outer diameter, n_tips)
    %   params   / struct / 1x1   Combustion parameters (from combustion_params())
    % OUTPUT
    %   lookup   / struct / 1x1   Geometry lookup with fields b, Ap, perim, b_max

    % Data extraction
    geometry_type = params.geometry.type;
    multiplier = params.geometry.multiplier;

    % Number of contour points
    if isfield(params.geometry, 'npoints') && ~isempty(params.geometry.npoints)
        npoints = params.geometry.npoints;
    elseif geometry_type == "star"
        npoints = (meshdata.n_tips + 1) + meshdata.n_tips * multiplier;
    else
        npoints = max(200, multiplier);
    end

    % Engine and time data
    ext_diameter = params.engine.ext_diameter;        % [m]
    tmax = params.time.tmax;                          % [s]
    time_output = 0:params.time.dt_output:tmax;       % [s]

    % MDF GEOMETRY LOOKUP 
    mdf_cfg = struct();
    mdf_cfg.mode = "lookup";
    mdf_cfg.geometry_type = geometry_type;
    mdf_cfg.meshdata = meshdata;
    mdf_cfg.npoints = npoints;
    mdf_cfg.casing_radius = ext_diameter / 2;
    mdf_cfg.t_vec = time_output;
    mdf_cfg.h = mdf_cfg.casing_radius / params.mdf.grid_divisions;
    mdf_cfg.use_bwdist = params.mdf.use_bwdist;
    mdf_cfg.use_parallel = params.mdf.use_parallel;
    mdf_cfg.do_plots = false;
    mdf_cfg.v_reg = 0;
    mdf_cfg.geom_opts.use_lookup = true;
    mdf_cfg.geom_opts.n_lookup = params.mdf.n_lookup;
    mdf_cfg.geom_opts.perimeter_from_area = params.mdf.perimeter_from_area;
    mdf_cfg.geom_opts.store_contours = false;
    mdf_cfg.geom_opts.smooth_prefix = params.mdf.smooth_prefix_frac * mdf_cfg.casing_radius;

    mdf = burnback_mdf_main(mdf_cfg);

    % Prepare the lookup arrays 
    lookup = clean_lookup(mdf);
end

%% Local functions

function lookup = clean_lookup(mdf)
    % clean_lookup
    % Clean the raw MDF output into a monotone, strictly increasing burnback
    % lookup and rebuild the port area from the burning perimeter.
    % INPUT
    %   mdf    / struct / 1x1   MDF solver output (b, port area, perimeter lookups)
    % OUTPUT
    %   lookup / struct / 1x1   Cleaned lookup with fields b, Ap, perim, b_max

    lookup.b = mdf.b_lookup(:);
    lookup.Ap = mdf.port_area_lookup(:);
    lookup.perim = mdf.perimeter_lookup(:);

    valid = isfinite(lookup.b) & isfinite(lookup.Ap) & isfinite(lookup.perim) & ...
        lookup.Ap > 0 & lookup.perim >= 0;
    lookup.b = lookup.b(valid);
    lookup.Ap = lookup.Ap(valid);
    lookup.perim = lookup.perim(valid);

    [lookup.b, unique_idx] = unique(lookup.b, 'stable');
    lookup.Ap = lookup.Ap(unique_idx);
    lookup.perim = lookup.perim(unique_idx);

    if numel(lookup.b) < 2
        error('MDF lookup must contain at least two valid burnback levels.');
    end
    if any(diff(lookup.b) <= 0)
        error('MDF burnback lookup must be strictly increasing.');
    end

    lookup.Ap_raw = lookup.Ap;
    lookup.Ap = lookup.Ap(1) + cumtrapz(lookup.b, lookup.perim);
    lookup.Ap = min(max(lookup.Ap, lookup.Ap_raw(1)), max(lookup.Ap_raw));
    lookup.dAp_db = lookup.perim;
    lookup.dperim_db = gradient(lookup.perim, lookup.b);
    lookup.b_max = lookup.b(end);
end
