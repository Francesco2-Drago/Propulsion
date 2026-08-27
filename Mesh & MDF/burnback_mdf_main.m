function out = burnback_mdf_main(cfg)
% BURNBACK_MDF_MAIN  Minimum Distance Function burnback driver.
%
%   out = burnback_mdf_main(cfg) runs the MDF burnback model using fields in
%   cfg. With no input it runs the original demo and plots the results.
%
%   cfg.mode = "lookup" returns geometry lookup data as a function of
%   cumulative burnback:
%       out.b_lookup
%       out.port_area_lookup
%       out.perimeter_lookup
%       out.half_contours_lookup

demo_mode = nargin < 1 || isempty(cfg);
if demo_mode
    close all;
    clc;
    cfg = struct();
end

cfg = default_burnback_main_cfg(cfg, demo_mode);

if isfield(cfg, 'half_mesh') && ~isempty(cfg.half_mesh)
    half_mesh = cfg.half_mesh;
else
    half_mesh = make_mesh0(cfg.geometry_type, cfg.meshdata, cfg.npoints, "cartesian");
end

% half_mesh is a (N x 2) array of [x,y] with y >= 0
max_r = max(sqrt(half_mesh(:,1).^2 + half_mesh(:,2).^2));
% Tolleranza per evitare falsi positivi dovuti a errori di arrotondamento
tol = 1e-10;
if cfg.casing_radius < max_r - tol
    error('Casing radius must be >= the outermost tip of the grain.');
elseif cfg.casing_radius < max_r
    % Se la differenza è solo numerica, correggiamo il raggio del casing
    cfg.casing_radius = max_r;
    warning('Casing radius adjusted by %.2e m to match grain tip.', max_r - cfg.casing_radius);
end

% fprintf('=== Grain Burnback MDF ===\n');
% t_sim = tic;
out = grain_burnback_mdf(half_mesh, cfg.v_reg, cfg.t_vec, cfg.casing_radius, cfg.h, ...
                         cfg.use_bwdist, cfg.use_parallel, cfg.geom_opts);
% fprintf('[SIM] Total time: %.2f s\n', toc(t_sim));

if cfg.mode == "lookup"
    return
end

if cfg.do_plots
    % fprintf('=== Visualizations ===\n');
    plot_phi0(out);
    plot_ballistics(out);
    animate_burnback(out, cfg.casing_radius);
end
end

function cfg = default_burnback_main_cfg(cfg, demo_mode)
if ~isfield(cfg, 'mode') || isempty(cfg.mode)
    cfg.mode = "history";
else
    cfg.mode = string(cfg.mode);
end

if ~isfield(cfg, 'geometry_type') || isempty(cfg.geometry_type)
    cfg.geometry_type = "star";
else
    cfg.geometry_type = string(cfg.geometry_type);
end

if ~isfield(cfg, 'meshdata') || isempty(cfg.meshdata)
    cfg.meshdata.inner_diameter = 0.7;   % [m]
    cfg.meshdata.outer_diameter = 1;     % [m]
    cfg.meshdata.n_tips = 8;
end

if ~isfield(cfg, 'multiplier') || isempty(cfg.multiplier)
    cfg.multiplier = 1000;
end

if ~isfield(cfg, 'npoints') || isempty(cfg.npoints)
    if cfg.geometry_type == "star"
        cfg.npoints = (cfg.meshdata.n_tips + 1) + cfg.meshdata.n_tips * cfg.multiplier;
    else
        cfg.npoints = max(200, cfg.multiplier);
    end
end

if ~isfield(cfg, 'casing_radius') || isempty(cfg.casing_radius)
    cfg.casing_radius = 1; % [m]
end

if ~isfield(cfg, 'v_reg') || isempty(cfg.v_reg)
    cfg.v_reg = 5e-3;   % [m/s] scalar, time-only function, or burn law struct
end

if ~isfield(cfg, 't_vec') || isempty(cfg.t_vec)
    cfg.t_vec = linspace(0, 300, 1001);
end

if ~isfield(cfg, 'h') || isempty(cfg.h)
    cfg.h = cfg.casing_radius / 1000;   % grid spacing [m]
end

if ~isfield(cfg, 'use_bwdist') || isempty(cfg.use_bwdist)
    cfg.use_bwdist = true;
end

if ~isfield(cfg, 'use_parallel') || isempty(cfg.use_parallel)
    cfg.use_parallel = false;
end

if ~isfield(cfg, 'geom_opts') || isempty(cfg.geom_opts)
    cfg.geom_opts = struct();
end
if ~isfield(cfg.geom_opts, 'use_lookup')
    cfg.geom_opts.use_lookup = true;
end
if ~isfield(cfg.geom_opts, 'n_lookup')
    cfg.geom_opts.n_lookup = numel(cfg.t_vec);
end
if ~isfield(cfg.geom_opts, 'perimeter_from_area')
    cfg.geom_opts.perimeter_from_area = false;
end
if ~isfield(cfg.geom_opts, 'store_grid')
    cfg.geom_opts.store_grid = cfg.mode ~= "lookup";
end
if ~isfield(cfg.geom_opts, 'store_contours')
    cfg.geom_opts.store_contours = true;
end
cfg.geom_opts.lookup_only = cfg.mode == "lookup";

if ~isfield(cfg, 'do_plots') || isempty(cfg.do_plots)
    cfg.do_plots = demo_mode;
end
end

%% ========================================================================
%  VISUALIZER FUNCTIONS
% =========================================================================

function plot_phi0(results)
% PLOT_PHI0  Colour map of phi0 on the HALF domain only.
% Ref: Zeriadtke et al. 2024, Fig.4-5.
% fprintf('[VIS] Plotting initial SDF field phi0(x,y) on half domain...\n');

phi0 = results.phi0;
X    = results.X;
Y    = results.Y;
xv   = results.xvec;
yv   = results.yvec;

figure('Name','Initial Signed Distance Field','Color','w');
contourf(X, Y, phi0, 50, 'LineColor','none');
hold on;
contour(X, Y, phi0, [0 0], 'k-', 'LineWidth', 2.0);
cb = colorbar;
cb.Label.String = '\phi_0 [m]';
colormap(gca, 'turbo');
xlabel('x [m]');
ylabel('y [m]');
title('\phi_0(x,y) - Initial Signed Distance Field (half domain)');
axis equal;
xlim([min(xv), max(xv)]);
ylim([min(yv), max(yv)]);
grid on;
end

function plot_ballistics(results)
% PLOT_BALLISTICS  Port area and burning perimeter vs time.
% fprintf('[VIS] Plotting internal ballistics geometry (area and perimeter vs time)...\n');

t  = results.t;
Ab = results.port_area;
P  = results.perimeter;

figure('Name','Internal Ballistics Geometry','Color','w', ...
       'Position',[100 100 680 440]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(t, Ab, 'b-', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('A_p [m^2]');
title('Port cross-sectional area');
grid on;
box on;

nexttile;
plot(t, P, 'r-', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('P_b [m]');
title('Burning perimeter');
grid on;
box on;
end

function animate_burnback(results, casing_radius, fps)
% ANIMATE_BURNBACK  Animated full-port evolution (half domain mirrored).
% fps defaults to n_frames/5 so animation lasts ~5 s.
if nargin < 3 || isempty(fps)
    fps = max(1, numel(results.t) / 5);
end
% fprintf('[VIS] Animating burnback evolution (fps = %.1f)...\n', fps);

t      = results.t;
contrs = results.half_contours;
n_t    = numel(t);
cmap   = parula(n_t);

th_c = linspace(0, 2*pi, 500)';
xc   = casing_radius * cos(th_c);
yc   = casing_radius * sin(th_c);

fig = figure('Name','Port Evolution','Color','w','Position',[150 150 640 640]);
ax  = axes('Parent', fig);
hold(ax,'on'); 
axis(ax,'equal'); 
grid(ax,'on'); 
box(ax,'on');
fill(ax, xc, yc, [0.22 0.22 0.22], 'EdgeColor','k', 'LineWidth',1.5);
xlim(ax, casing_radius * [-1.12, 1.12]);
ylim(ax, casing_radius * [-1.12, 1.12]);
xlabel(ax,'x [m]');
ylabel(ax,'y [m]');
title(ax,'Grain port evolution');
colormap(ax, parula);
clim(ax, [t(1), t(end)]);
cb = colorbar(ax); 
cb.Label.String = 'Time [s]';

frame_h = gobjects(0);
for k = 1:n_t
    delete(frame_h);
    polys = contrs{k};
    frame_h = gobjects(numel(polys) + 1, 1);
    col   = cmap(k,:);
    for p = 1:numel(polys)
        xy_top  = polys{p};
        xy_bot  = [xy_top(:,1), -xy_top(:,2)];
        xy_full = [xy_top; xy_bot(end:-1:1,:)];
        h_p = fill(ax, xy_full(:,1), xy_full(:,2), col, ...
                   'FaceAlpha',0.88, 'EdgeColor',col*0.7, 'LineWidth',0.8);
        frame_h(p) = h_p;
    end
    h_txt = text(ax, casing_radius*0.62, casing_radius*0.95, ...
                 sprintf('t = %.2f s', t(k)), ...
                 'FontSize',10, 'Color','w', 'HorizontalAlignment','right');
    frame_h(end) = h_txt;
    drawnow;
    pause(1/fps);
end
end

%% ========================================================================
%  GRAIN_BURNBACK_MDF
% =========================================================================
function results = grain_burnback_mdf(half_mesh_xy, v_reg, t_vec, ...
                                      casing_radius, h, use_bwdist, use_parallel, geom_opts)
% GRAIN_BURNBACK_MDF  2-D half-domain burnback via Minimum Distance Function.
%
% References:
%   [1] Willcox et al., J. Propulsion and Power, 23(2), 2007.
%   [2] Zeriadtke et al., Aerospace, 11(2), 103, 2024.
%
%   For uniform regression, phi_t = phi0 - d_b and phi_t < 0 <=> phi0 < d_b,
%   instead of calculating phi_t = phi0 - d_b as a full matrix each step,
%   pass phi0 directly to contourc with the shifted iso-level d_b:
%       contourc(phi0_sub, [d_b d_b])
%   If v_reg is coupled to G_ox but uniform over the section, only the scalar
%   cumulative burnback d_b must be advanced serially. Area and perimeter are
%   read from a precomputed d_b lookup table.
%   A spatial v_reg(X,Y,t) still requires direct level-set evolution.

if nargin < 6
    use_bwdist  = true;
end
if nargin < 7
    use_parallel = false;
end
if nargin < 8 || isempty(geom_opts)
    geom_opts = struct();
end
geom_opts = default_geom_opts(geom_opts);

v_is_scalar       = isnumeric(v_reg) && isscalar(v_reg);
v_is_rate_history = isnumeric(v_reg) && isvector(v_reg) && numel(v_reg) == numel(t_vec) && ~v_is_scalar;
v_is_func         = isa(v_reg, 'function_handle');
v_is_burn_law     = is_burn_law(v_reg);
v_nargin          = function_handle_nargin(v_reg);
v_is_uniform_func = v_is_func && (v_nargin == 1 || v_nargin == 4 || v_nargin < 0);
v_is_spatial_func = v_is_func && v_nargin == 3;

% parfor is only valid for independent geometry levels, not for the coupled
% time march where G_ox depends on the previous port area.
if use_parallel && ~(v_is_scalar || v_is_rate_history || v_is_uniform_func || v_is_burn_law)
    warning('use_parallel = true requires independent geometry levels. Switching to serial.');
    use_parallel = false;
end

% =========================================================================
% 1. GRID
% =========================================================================
pad  = max(casing_radius * 0.05, 5*h);
xmin = min(min(half_mesh_xy(:,1)) - pad, -casing_radius - pad);
xmax = max(max(half_mesh_xy(:,1)) + pad, casing_radius + pad);
ymin = 0;
ymax = max(max(half_mesh_xy(:,2)) + pad, casing_radius + pad);

xvec = (xmin:h:xmax)';
yvec = (ymin:h:ymax)';
Nx = numel(xvec);
Ny = numel(yvec);
need_grid_coordinates = geom_opts.store_grid || ~use_bwdist || v_is_spatial_func;
if need_grid_coordinates
    [X, Y] = meshgrid(xvec, yvec);
else
    X = [];
    Y = [];
end
% fprintf('[MDF] Grid (half): %d x %d = %d points\n', Ny, Nx, Ny*Nx);

% =========================================================================
% 2. INITIAL SDF phi0
%    phi0 < 0 => gas port interior
%    phi0 > 0 => solid / outside                         [Ref 1, Sec.II]
% =========================================================================
% fprintf('[MDF] Building phi0 ');
% t_sdf = tic;

hm        = half_mesh_xy;
closed_xy = [hm; hm(end,1), 0; hm(1,1), 0];
mask_port = scanline_fill(closed_xy(:,1), closed_xy(:,2), xvec, yvec, Nx, Ny);

if use_bwdist
    % fprintf('(bwdist)...\n');
    dist_out = double(bwdist( mask_port)) * h;
    dist_in = double(bwdist(~mask_port)) * h;
    phi0 = dist_out - dist_in;
else
    % fprintf('(custom SDF)...\n');
    phi0 = sdf_poly(X, Y, hm);
    phi0(mask_port) = -phi0(mask_port);
end

% fprintf('[MDF] phi0 done in %.2f s\n', toc(t_sdf));

% =========================================================================
% 3. CASING MASK & BURNOUT BOUNDARY
%    Outside-casing phi values are forced positive (never burn).
%    [Ref 1, Sec.III.C.2]
% =========================================================================
inside_casing  = (yvec.^2 + (xvec').^2) <= casing_radius^2;
outside_casing = ~inside_casing;
phi0(outside_casing) = abs(phi0(outside_casing));   % enforce wall

half_casing_area = pi * casing_radius^2 / 2;
h2               = h * h;
R                = casing_radius;
R2               = R * R;
burn_limit2      = (R - 0.5*h)^2;
tol_P            = max(h, 1e-12);

% Pre-index inside-casing cells for fast area count
ic_idx   = find(inside_casing);

% Row/col masks for bounding-box computation (reused every step)
xvec_row = xvec';
yvec_row = yvec';

% Bounding boxes for contour extraction can be found from row/column minima
% of phi0 inside the casing. This avoids rebuilding active = phi0 < b for
% every contour level.
phi0_inside = phi0;
phi0_inside(~inside_casing) = inf;
row_min_phi = min(phi0_inside, [], 2);
col_min_phi = min(phi0_inside, [], 1);
clear phi0_inside

% =========================================================================
% 4. SURFACE EVOLUTION
%
%  Uniform combustion:
%    phi_t = phi0 - b(t), so the geometry is a single-variable function of
%    cumulative burnback b. If G_ox changes, only b(t) is marched serially;
%    area/perimeter are interpolated from a precomputed b lookup table.
%
%  Spatially non-uniform combustion:
%    db is a matrix, phi_t must be rebuilt, and no scalar geometry lookup is
%    valid.
%                                                           [Ref 1, Sec.III]
% =========================================================================
n_t           = numel(t_vec);
port_area_h   = nan(n_t, 1);
perimeter_h   = nan(n_t, 1);
burnback_h    = nan(n_t, 1);
reg_rate_h    = nan(n_t, 1);
gox_h         = nan(n_t, 1);
half_contours = cell(n_t, 1);

% Pre-extract phi0 values at inside-casing indices.
phi0_ic        = phi0(ic_idx);
max_burnback_h = max(phi0_ic) + h;

if isfield(geom_opts, 'lookup_only') && geom_opts.lookup_only
    lookup_b = linspace(0, max_burnback_h, geom_opts.n_lookup).';
    lookup = build_geometry_lookup(lookup_b, phi0, phi0_ic, ...
        xvec_row, yvec_row, row_min_phi, col_min_phi, h2, half_casing_area, R, R2, burn_limit2, ...
        max_burnback_h, use_parallel, geom_opts, hm);
    % fprintf('[MDF] Geometry lookup: %d burnback levels\n', numel(lookup_b));

    results.t                     = t_vec(:)';
    results.port_area             = [];
    results.perimeter             = [];
    results.burnback              = [];
    results.reg_rate              = [];
    results.gox                   = [];
    results.half_contours         = {};
    results.phi0                  = phi0;
    results.X                     = X;
    results.Y                     = Y;
    results.xvec                  = xvec;
    results.yvec                  = yvec;
    results.b_lookup              = lookup.burnback(:);
    results.port_area_lookup      = 2 * lookup.port_area_h(:);
    results.perimeter_lookup      = 2 * lookup.perimeter_h(:);
    results.half_contours_lookup  = lookup.half_contours(:);
    results.max_burnback          = max_burnback_h;
    return
end

if v_is_scalar || v_is_rate_history
    % ---- KNOWN UNIFORM BURNBACK HISTORY ------------------------------
    if v_is_scalar
        reg_rate_h(:) = v_reg;
        burnback_h = v_reg * (t_vec(:) - t_vec(1));
    else
        reg_rate_h = v_reg(:);
        burnback_h(1) = 0;
        for k = 1:n_t-1
            burnback_h(k+1) = burnback_h(k) + reg_rate_h(k) * (t_vec(k+1) - t_vec(k));
        end
    end

    [port_area_h, perimeter_h, half_contours] = build_geometry_history( ...
        burnback_h, phi0, phi0_ic, xvec_row, yvec_row, row_min_phi, col_min_phi, ...
        h2, half_casing_area, R, R2, burn_limit2, max_burnback_h, ...
        use_parallel);

    ib = find(perimeter_h <= tol_P, 1);
    if ~isempty(ib)
        t_vec         = t_vec(1:ib);
        port_area_h   = port_area_h(1:ib);
        perimeter_h   = perimeter_h(1:ib);
        burnback_h    = burnback_h(1:ib);
        reg_rate_h    = reg_rate_h(1:ib);
        gox_h         = gox_h(1:ib);
        half_contours = half_contours(1:ib);
        % fprintf('[MDF] Burnout at t = %.4f\n', t_vec(end));
    end

elseif v_is_uniform_func || v_is_burn_law
    % ---- Gox-COUPLED UNIFORM REGRESSION ------------------------------
    if geom_opts.use_lookup
        lookup_b = linspace(0, max_burnback_h, geom_opts.n_lookup).';
        lookup = build_geometry_lookup(lookup_b, phi0, phi0_ic, ...
            xvec_row, yvec_row, row_min_phi, col_min_phi, h2, half_casing_area, R, R2, burn_limit2, ...
            max_burnback_h, use_parallel, geom_opts, hm);
        % fprintf('[MDF] Geometry lookup: %d burnback levels\n', numel(lookup_b));
    else
        lookup = [];
    end

    b_now = 0;
    for k = 1:n_t
        if geom_opts.use_lookup
            [A_h, P_h, polys] = interpolate_geometry_lookup(lookup, b_now);
        else
            [A_h, P_h, polys] = geometry_at_level(b_now, phi0, phi0_ic, ...
                xvec_row, yvec_row, row_min_phi, col_min_phi, h2, half_casing_area, ...
                R, R2, burn_limit2, max_burnback_h);
        end

        port_area_h(k)   = A_h;
        perimeter_h(k)   = P_h;
        burnback_h(k)    = b_now;
        half_contours{k} = polys;

        if P_h <= tol_P
            t_vec         = t_vec(1:k);
            port_area_h   = port_area_h(1:k);
            perimeter_h   = perimeter_h(1:k);
            burnback_h    = burnback_h(1:k);
            reg_rate_h    = reg_rate_h(1:k);
            gox_h         = gox_h(1:k);
            half_contours = half_contours(1:k);
            % fprintf('[MDF] Burnout at t = %.4f\n', t_vec(k));
            break
        end

        [reg_rate_h(k), gox_h(k)] = eval_uniform_rate(v_reg, t_vec(k), 2*A_h, 2*P_h, b_now);
        if k < n_t
            b_now = b_now + reg_rate_h(k) * (t_vec(k+1) - t_vec(k));
        end
    end

elseif v_is_spatial_func || isnumeric(v_reg)
    % ---- SPATIALLY NON-UNIFORM REGRESSION ----------------------------
    db = zeros(Ny, Nx);

    for k = 1:n_t
        dt = t_vec(k) - (k > 1) * t_vec(k-1);
        if v_is_spatial_func
            db = db + v_reg(X, Y, t_vec(k)) * dt;
        else
            db = db + v_reg * dt;
        end
        phi_t = phi0 - db;
        % Re-enforce boundary (db may have grown non-uniformly)
        phi_t(outside_casing) = max(phi_t(outside_casing), 1e-10);

        n_gas = sum(phi_t(ic_idx) < 0);
        A_h   = min(h2 * n_gas, half_casing_area);
        port_area_h(k) = A_h;
        burnback_h(k)  = mean(db(ic_idx));

        if n_gas == 0
            half_contours{k} = {};
            perimeter_h(k)   = 0;
            break
        end

        active = phi_t < 0;
        rows   = any(active, 2);
        cols   = any(active, 1);
        rmin   = max(1,  find(rows,1,'first') - 1);
        rmax   = min(Ny, find(rows,1,'last')  + 1);
        cmin   = max(1,  find(cols,1,'first') - 1);
        cmax   = min(Nx, find(cols,1,'last')  + 1);

        C = contourc(xvec_row(cmin:cmax), yvec_row(rmin:rmax), ...
                     phi_t(rmin:rmax, cmin:cmax), [0 0]);

        [polys, P_h] = extract_contour(C, R, R2, burn_limit2);
        half_contours{k} = polys;
        perimeter_h(k)   = P_h;

        if P_h <= tol_P
            t_vec         = t_vec(1:k);
            port_area_h   = port_area_h(1:k);
            perimeter_h   = perimeter_h(1:k);
            burnback_h    = burnback_h(1:k);
            reg_rate_h    = reg_rate_h(1:k);
            gox_h         = gox_h(1:k);
            half_contours = half_contours(1:k);
            % fprintf('[MDF] Burnout at t = %.4f\n', t_vec(k));
            break
        end
    end
else
    error('Unsupported v_reg type.');
end

results.t             = t_vec(:)';
results.port_area     = 2 * port_area_h;
results.perimeter     = 2 * perimeter_h;
results.burnback      = burnback_h(:)';
results.reg_rate      = reg_rate_h(:)';
results.gox           = gox_h(:)';
results.half_contours = half_contours;
results.phi0          = phi0;
results.X             = X;
results.Y             = Y;
results.xvec          = xvec;
results.yvec          = yvec;
results.max_burnback  = max_burnback_h;
if exist('lookup', 'var') && ~isempty(lookup)
    results.b_lookup             = lookup.burnback(:);
    results.port_area_lookup     = 2 * lookup.port_area_h(:);
    results.perimeter_lookup     = 2 * lookup.perimeter_h(:);
    results.half_contours_lookup = lookup.half_contours(:);
end
end


%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function geom_opts = default_geom_opts(geom_opts)
if ~isfield(geom_opts, 'use_lookup')
    geom_opts.use_lookup = true;
end
if ~isfield(geom_opts, 'n_lookup')
    geom_opts.n_lookup = 1001;
end
if ~isfield(geom_opts, 'store_grid')
    geom_opts.store_grid = true;
end
if ~isfield(geom_opts, 'store_contours')
    geom_opts.store_contours = true;
end
geom_opts.n_lookup = max(2, round(geom_opts.n_lookup));
end

function n = function_handle_nargin(fh)
if ~isa(fh, 'function_handle')
    n = NaN;
    return
end
try
    n = nargin(fh);
catch
    n = -1;
end
end

function tf = is_burn_law(rate_model)
tf = isstruct(rate_model) && ...
     all(isfield(rate_model, {'mdot_ox', 'a_rf', 'n_rf'}));
end

function [v, gox] = eval_uniform_rate(rate_model, t, Ap, Pb, b)
gox = NaN;
if is_burn_law(rate_model)
    mdot_ox = eval_time_value(rate_model.mdot_ox, t);
    gox = mdot_ox / Ap;
    v = rate_model.a_rf * gox^rate_model.n_rf;
else
    n = function_handle_nargin(rate_model);
    if n == 1
        v = rate_model(t);
    else
        v = rate_model(t, Ap, Pb, b);
    end
end
if ~isscalar(v) || ~isfinite(v)
    error('Uniform v_reg function must return one finite scalar value.');
end
if ~isscalar(gox) || ~(isfinite(gox) || isnan(gox))
    error('G_ox must be one finite scalar value.');
end
end

function value = eval_time_value(value_or_fun, t)
if isa(value_or_fun, 'function_handle')
    value = value_or_fun(t);
else
    value = value_or_fun;
end
if ~isscalar(value) || ~isfinite(value)
    error('Time-dependent burn-law input must return one finite scalar value.');
end
end

function lookup = build_geometry_lookup(db_levels, phi0, phi0_ic, ...
                                        xvec_row, yvec_row, row_min_phi, col_min_phi, ...
                                        h2, half_casing_area, R, R2, burn_limit2, ...
                                        max_burnback_h, use_parallel, geom_opts, initial_contour)
if nargin < 15 || isempty(geom_opts)
    geom_opts = struct();
end
if nargin < 16
    initial_contour = [];
end
if isfield(geom_opts, 'perimeter_from_area') && geom_opts.perimeter_from_area
    db_levels = db_levels(:);
    port_area_h = area_from_levels(phi0_ic, db_levels, h2, half_casing_area);
    port_area_h = port_area_h(:);
    perimeter_h = max(gradient(port_area_h, db_levels), 0);
    half_contours = cell(numel(db_levels), 1);
    if isempty(initial_contour)
        [half_contours{1}, P0_h] = contour_at_level(db_levels(1), phi0, ...
            xvec_row, yvec_row, row_min_phi, col_min_phi, R, R2, burn_limit2, max_burnback_h);
    else
        half_contours{1} = {initial_contour};
        P0_h = sum(hypot(diff(initial_contour(:,1)), diff(initial_contour(:,2))));
    end
    if isfinite(P0_h) && P0_h > 0
        perimeter_h(1) = P0_h;
    end
    positive = find(isfinite(perimeter_h) & perimeter_h > 0);
    if ~isempty(positive)
        first_positive = positive(1);
        last_positive = positive(end);
        perimeter_h(1:first_positive-1) = perimeter_h(first_positive);
        if last_positive > first_positive
            idx = first_positive:last_positive;
            perimeter_h(idx) = interp1(db_levels(positive), perimeter_h(positive), ...
                db_levels(idx), 'linear');
        end
    end
else
    store_contours = geom_opts.store_contours;
    [port_area_h, perimeter_h, half_contours] = build_geometry_history( ...
        db_levels, phi0, phi0_ic, xvec_row, yvec_row, row_min_phi, col_min_phi, ...
        h2, half_casing_area, R, R2, burn_limit2, max_burnback_h, use_parallel, ...
        store_contours);
end

lookup.burnback      = db_levels(:);
lookup.port_area_h   = port_area_h(:);
lookup.perimeter_h   = perimeter_h(:);
lookup.half_contours = half_contours(:);

% LUT smoothing: replace geometry values in [b=0, b=smooth_prefix]
% with interpolated values anchored at the exact initial geometry and 
% the original LUT value at the prefix end
if isfield(geom_opts, 'smooth_prefix') && isscalar(geom_opts.smooth_prefix) && geom_opts.smooth_prefix > 0
    sp = geom_opts.smooth_prefix;
    b = db_levels(:);
    % only apply if prefix extends beyond the first sample
    if sp > b(1)
        % last index inside the prefix (if none, take first)
        idx_last = find(b <= sp, 1, 'last');
        if isempty(idx_last)
            idx_last = 1;
        end

        % try to use exact initial geometry from mesh (initial_contour)
        if exist('initial_contour', 'var') && ~isempty(initial_contour)
            try
                [P0_full, A0_full] = mesh_eval(initial_contour, 'cartesian');
                A0_h = A0_full / 2; % convert to half-domain value used internally
                P0_h = P0_full / 2;
            catch
                % fall back to the LUT first value
                A0_h = port_area_h(1);
                P0_h = perimeter_h(1);
            end
        else
            A0_h = port_area_h(1);
            P0_h = perimeter_h(1);
        end

        % anchors: exact initial and the LUT value at the prefix end.
        anchors_b = [b(1); b(idx_last)];
        anchors_A = [A0_h; port_area_h(idx_last)];
        anchors_P = [P0_h; perimeter_h(idx_last)];

        % shape-preserving interpolant (pchip) on the prefix nodes
        new_b = b(1:idx_last);
        newA = interp1(anchors_b, anchors_A, new_b, 'pchip');
        newP = interp1(anchors_b, anchors_P, new_b, 'pchip');

        % port area must be non-decreasing with burnback
        % enforce monotonicity to avoid tiny decreases for interpolation/numerical issues
        newA = cummax(newA);

        % replace prefix values in the LUT and update lookup outputs
        port_area_h(1:idx_last) = newA;
        perimeter_h(1:idx_last) = newP;
        lookup.port_area_h = port_area_h(:);
        lookup.perimeter_h = perimeter_h(:);
    end
end
end

function [A_h, P_h, polys] = interpolate_geometry_lookup(lookup, b)
b = min(max(b, lookup.burnback(1)), lookup.burnback(end));
A_h = interp1(lookup.burnback, lookup.port_area_h, b, 'linear');
P_h = interp1(lookup.burnback, lookup.perimeter_h, b, 'linear');
P_h = max(P_h, 0);

[~, idx] = min(abs(lookup.burnback - b));
polys = lookup.half_contours{idx};
end

function [port_area_h, perimeter_h, half_contours] = build_geometry_history( ...
    db_levels, phi0, phi0_ic, xvec_row, yvec_row, row_min_phi, col_min_phi, ...
    h2, half_casing_area, R, R2, burn_limit2, max_burnback_h, use_parallel, store_contours)

if nargin < 15
    store_contours = true;
end

db_levels     = db_levels(:);
n_levels      = numel(db_levels);
port_area_h   = area_from_levels(phi0_ic, db_levels, h2, half_casing_area);
port_area_h   = port_area_h(:);
perimeter_h   = nan(n_levels, 1);
half_contours = cell(n_levels, 1);

if use_parallel
    parfor k = 1:n_levels
        [polys, P_h] = contour_at_level(db_levels(k), phi0, ...
            xvec_row, yvec_row, row_min_phi, col_min_phi, R, R2, burn_limit2, max_burnback_h);
        perimeter_h(k)   = P_h;
        if store_contours || k == 1
            half_contours{k} = polys;
        else
            half_contours{k} = {};
        end
    end
else
    for k = 1:n_levels
        [polys, P_h] = contour_at_level(db_levels(k), phi0, ...
            xvec_row, yvec_row, row_min_phi, col_min_phi, R, R2, burn_limit2, max_burnback_h);
        perimeter_h(k)   = P_h;
        if store_contours || k == 1
            half_contours{k} = polys;
        else
            half_contours{k} = {};
        end
    end
end
end

function area_h = area_from_levels(phi0_ic, db_levels, h2, half_casing_area)
out_size = size(db_levels);
db_col = db_levels(:);
[db_sorted, sort_idx] = sort(db_col);
[db_unique, ~, unique_idx] = unique(db_sorted);
edges = [-inf; db_unique; inf];
bin_counts = histcounts(phi0_ic, edges);
area_unique = min(h2 * cumsum(bin_counts(1:end-1)).', half_casing_area);
area_sorted = area_unique(unique_idx);
area_col = nan(size(db_col));
area_col(sort_idx) = area_sorted;
area_h = reshape(area_col, out_size);
end

function [A_h, P_h, polys] = geometry_at_level(b, phi0, phi0_ic, ...
                                               xvec_row, yvec_row, row_min_phi, col_min_phi, ...
                                               h2, half_casing_area, ...
                                               R, R2, burn_limit2, max_burnback_h)
A_h = area_from_levels(phi0_ic, b, h2, half_casing_area);
[polys, P_h] = contour_at_level(b, phi0, xvec_row, yvec_row, ...
    row_min_phi, col_min_phi, R, R2, burn_limit2, max_burnback_h);
end

function [polys, P_h] = contour_at_level(b, phi0, xvec_row, yvec_row, ...
                                         row_min_phi, col_min_phi, R, R2, burn_limit2, max_burnback_h)
P_h = 0;
polys = {};

if b >= max_burnback_h
    return
end

rows = row_min_phi < b;
cols = col_min_phi < b;
if ~any(rows) || ~any(cols)
    return
end

rmin = max(1,  find(rows, 1, 'first') - 1);
rmax = min(size(phi0, 1), find(rows, 1, 'last')  + 1);
cmin = max(1,  find(cols, 1, 'first') - 1);
cmax = min(size(phi0, 2), find(cols, 1, 'last')  + 1);

C = contourc(xvec_row(cmin:cmax), yvec_row(rmin:rmax), ...
             phi0(rmin:rmax, cmin:cmax), [b b]);
[polys, P_h] = extract_contour(C, R, R2, burn_limit2);
end

function [polys, P_h] = extract_contour(C, R, R2, burn_limit2)
% extract contours, clamp to casing, accumulate perimeter
P_h    = 0;
polys  = cell(1, size(C, 2));
np_out = 0;
idx    = 1;
n_C    = size(C, 2);
while idx < n_C
    np  = C(2, idx);
    pts = C(:, idx+1 : idx+np)';

    % Clamp escaped vertices back onto casing  [Ref 1, Sec.III.C.2]
    r2_pts  = pts(:,1).^2 + pts(:,2).^2;
    outside = r2_pts > R2;
    if any(outside)
        sc = R ./ sqrt(r2_pts(outside));
        pts(outside,1) = pts(outside,1) .* sc;
        pts(outside,2) = pts(outside,2) .* sc;
    end

    np_out = np_out + 1;
    polys{np_out} = pts;

    % Perimeter: only segments whose midpoint is strictly inside casing
    dx  = diff(pts(:,1));
    dy  = diff(pts(:,2));
    xm  = pts(1:end-1,1) + 0.5*dx;
    ym  = pts(1:end-1,2) + 0.5*dy;
    rm2 = xm.^2 + ym.^2;
    burning = rm2 <= burn_limit2;
    P_h = P_h + sum(hypot(dx(burning), dy(burning)));

    idx = idx + np + 1;
end
polys = polys(1:np_out);
end

function mask = scanline_fill(px, py, xvec, yvec, Nx, Ny)
% scanline even-odd fill
mask = false(Ny, Nx);
n_v  = numel(px);
xmin = xvec(1);
h    = xvec(2) - xvec(1);
for row = 1:Ny
    y_row     = yvec(row);
    crossings = zeros(1, n_v);
    nc = 0;
    for e = 1:n_v-1
        y1 = py(e);  y2 = py(e+1);
        if (y1 <= y_row && y2 > y_row) || (y2 <= y_row && y1 > y_row)
            nc = nc + 1;
            crossings(nc) = px(e) + (y_row - y1)/(y2 - y1) * (px(e+1) - px(e));
        end
    end
    if nc == 0, continue; end
    crossings = sort(crossings(1:nc));
    for c = 1:2:nc-1
        j1 = max(1,  ceil( (crossings(c)   - xmin)/h) + 1);
        j2 = min(Nx, floor((crossings(c+1) - xmin)/h) + 1);
        if j1 <= j2
            mask(row, j1:j2) = true;
        end
    end
end
end

function d = sdf_poly(X, Y, poly_xy)
% custom SDF (no toolbox) (slow as fuck)
Xf    = X(:);  Yf = Y(:);
n_seg = size(poly_xy,1) - 1;
d2min = inf(numel(Xf), 1);
for s = 1:n_seg
    ax = poly_xy(s,1);    ay = poly_xy(s,2);
    bx = poly_xy(s+1,1);  by = poly_xy(s+1,2);
    abx = bx-ax;  aby = by-ay;
    ab2 = abx^2 + aby^2;
    if ab2 < 1e-24, continue; end
    t     = max(0, min(1, ((Xf-ax)*abx + (Yf-ay)*aby) / ab2));
    cx    = ax + t*abx;  cy = ay + t*aby;
    d2min = min(d2min, (Xf-cx).^2 + (Yf-cy).^2);
end
d = reshape(sqrt(d2min), size(X));
end
