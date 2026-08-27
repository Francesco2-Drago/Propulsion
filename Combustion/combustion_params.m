function params = combustion_params()
%COMBUSTION_PARAMS Baseline inputs for a near-neutral 50 kN burn.
%
% This case is intended as a practical starting point, not as a final
% optimum. The shallow star behaves close to a circular port, so the burning
% perimeter changes mildly during most of the burn and tails off near the
% case at the end.

%% Mission target
% Constraint values C6, C7, C8: read from optimizationConstraints, which is the
% single source of truth for every threshold. Do not write numbers here.
K_constraints = constraint_table();
params.target.thrust = K_constraints.C7.lo;       % [N]   C7
params.target.p_chamber = K_constraints.C8.lo;    % [Pa]  C8
params.target.burn_time = K_constraints.C6.lo;    % [s]   C6

%% Geometry

params.geometry.type = "cylinder";
params.geometry.meshdata.diameter = 0.194;         % [m]
params.geometry.multiplier = 150;

%% Fuel: HTPB
params.fuel.rho_f = 920;       % [kg/m3]
params.fuel.a_rf = 0.027e-3;   % [m/s] / [kg/(m2 s)]^n
params.fuel.n_rf = 0.75;       % [-]

%% Combustion and thermochemistry
% The oxidizer is a FREE parameter (C14 withdrawn: "GOX" in the assignment is
% the oxidizer flux G_ox, not gaseous oxygen). This is only the default for a
% standalone run; main_optimization enumerates the candidates. Gaseous O2 is
% excluded on storage grounds, see C14 in optimizationConstraints.
params.combustion.oxidizer_lookup = "H2O2_90";
params.combustion.mdot_ox = 12.435;   % [kg/s]
params.combustion.of_min = K_constraints.C5.lo;   % [-]  C5, CEA domain lower bound
params.combustion.of_max = K_constraints.C5.hi;   % [-]  C5, CEA domain upper bound
% Chamber-pressure extent of the CEA tables. These are ODE event bounds, i.e.
% the validity range of the lookup, NOT a design constraint: the design limit
% on pressure is C8.
params.combustion.p_min = 101325;     % [Pa]
params.combustion.p_max = 99e5;       % [Pa]

%% Engine
params.engine.ext_diameter = 0.654;       % [m]
params.engine.chamber_length = 2.548;     % [m]
params.engine.throat_diameter = 0.1301;   % [m]
params.engine.eps = 200;                  % [-]
params.engine.pamb = 0;                   % [Pa]

%% Time integration
params.time.tmax = params.target.burn_time;    % [s]
params.time.dt_output = 0.25;                  % [s]
params.time.fine_ode = true;

%% MDF lookup
params.mdf.grid_divisions = 700;        % h = casing_radius / grid_divisions
params.mdf.n_lookup = 601;              % burnback levels
params.mdf.use_bwdist = true;
params.mdf.use_parallel = false;
params.mdf.perimeter_from_area = false;
params.mdf.smooth_prefix_frac = 0.015;  % [-] initial LUT smoothing, as a FRACTION
                                        %     of the casing radius (0 = disabled).
                                        %     Must be relative: an absolute length
                                        %     would break the scale invariance the
                                        %     normalized shape phase relies on.

%% Output
params.output.do_plots = usejava("desktop");
params.output.do_animation = params.output.do_plots;
params.output.animation_max_frames = 120;
params.output.animation_pause = 0.03;
params.output.animation_autoplay = false;

end

%% FUNCTIONS

function K = constraint_table()
    % constraint_table
    % Fetch the id-keyed constraint table, making sure the repository is on the
    % path first: this file must keep working when it is called on its own.
    % INPUT
    %   None
    % OUTPUT
    %   K / struct / 1x1   Constraints keyed by id (see optimizationConstraints)

    if exist("optimizationConstraints", "file") ~= 2
        repoRoot = fileparts(fileparts(mfilename("fullpath")));
        addpath(genpath(repoRoot));
    end
    [~, K] = optimizationConstraints();
end