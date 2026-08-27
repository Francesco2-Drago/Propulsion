clear; clc; close all;

% Checkpoint Tappa 1: verifica build_shape_lookup + confronto perimeter_from_area

%% Parameters and configurations
repoRoot = fileparts(pwd);           % pwd = ...\Optimization
addpath(genpath(repoRoot));

params = combustion_params();
params.geometry.type = "star";

% Casing radius and star shape defined as ratios of R_c (stay < 1 to fit)
R_c = params.engine.ext_diameter / 2;   % [m] casing radius

ri_ratio = 0.72;    % inner apothem / R_c        [-]
re_ratio = 0.95;    % tip radius / R_c   (< 1)    [-]
N        = 10;      % number of tips              [-]

meshdata.inner_diameter = 2 * ri_ratio * R_c;    % [m]
meshdata.outer_diameter = 2 * re_ratio * R_c;    % [m]
meshdata.n_tips = N;

% Containment check (same conditions used by the cost function pre-check)
r_e = re_ratio * R_c;                        % [m] tip radius
r_vertex = (ri_ratio * R_c) / cos(pi / N);   % [m] inner polygon vertex
assert(r_e <= R_c, "Tip radius %.3f m exceeds casing %.3f m", r_e, R_c);
assert(r_vertex <= R_c, "Vertex radius %.3f m exceeds casing %.3f m", r_vertex, R_c);
assert(r_e >= r_vertex, "Invalid star: tip radius below vertex radius");

%% Lookup build - contour perimeter vs perimeter-from-area
% Same shape, two ways of computing the burning perimeter

params_contour = params;
params_contour.mdf.perimeter_from_area = false;
lookup_contour = build_shape_lookup(meshdata, params_contour);

params_grad = params;
params_grad.mdf.perimeter_from_area = true;
lookup_grad = build_shape_lookup(meshdata, params_grad);

% Console summary
fprintf("Casing radius R_c = %.4f m\n", R_c);
fprintf("--- perimeter_from_area = false ---\n");
fprintf("\tAp(0)    = %.6e m^2\n", lookup_contour.Ap(1));
fprintf("\tperim(0) = %.6e m\n",   lookup_contour.perim(1));
fprintf("\tb_max    = %.6e m\n",   lookup_contour.b_max);
fprintf("--- perimeter_from_area = true  ---\n");
fprintf("\tAp(0)    = %.6e m^2\n", lookup_grad.Ap(1));
fprintf("\tperim(0) = %.6e m\n",   lookup_grad.perim(1));
fprintf("\tb_max    = %.6e m\n",   lookup_grad.b_max);

%% Plot comparison
figure('Name', 'Shape lookup check', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(lookup_contour.b, lookup_contour.Ap, 'LineWidth', 1.5); hold on
plot(lookup_grad.b, lookup_grad.Ap, '--', 'LineWidth', 1.5);
grid on
xlabel('b [m]'); ylabel('A_p [m^2]');
legend("contour", "from area", "Location", "best");
title('Port area');

nexttile
plot(lookup_contour.b, lookup_contour.perim, 'LineWidth', 1.5); hold on
plot(lookup_grad.b, lookup_grad.perim, '--', 'LineWidth', 1.5);
grid on
xlabel('b [m]'); ylabel('P_b [m]');
legend("contour", "from area", "Location", "best");
title('Burning perimeter');
