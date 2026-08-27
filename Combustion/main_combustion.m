%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LABORATORIO DI PROPULSIONE AEROSPAZIALE                %
%                              A.A. 2025/2026                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CODE
% Hybrid rocket engine simulation: build the grain-geometry lookup (port area
% and burning perimeter vs burnback) with the Marker-Distance-Field solver,
% find the steady-state chamber pressure, integrate the coupled pressure /
% burnback ODE over the burn, and post-process the engine performance
% (thrust, specific impulse, mixture ratio, ...).
    
close all
clear
clc

% Make every sub-folder of the repository visible to the script.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

%% STEP 0: data
% Load the parameter struct and unpack geometry, fuel, combustion and engine inputs.

tic
params = combustion_params();

% ------------ GEOMETRY DATA ------------
% Grain geometry descriptors and number of points used to discretise it.
geometry_type = params.geometry.type;
meshdata = params.geometry.meshdata;
multiplier = params.geometry.multiplier;
if isfield(params.geometry, 'npoints') && ~isempty(params.geometry.npoints)
    npoints = params.geometry.npoints;
elseif geometry_type == "star"
    npoints = (meshdata.n_tips + 1) + meshdata.n_tips * multiplier;
else
    npoints = max(200, multiplier);
end

% ------------ FUEL DATA ------------
% FUEL = HTPB
% regression rate rf = a*Gox^n
rho_f = params.fuel.rho_f;        % Fuel density [kg/m3]
a_rf = params.fuel.a_rf;          % Regression law coefficient [m/s]
n_rf = params.fuel.n_rf;          % Regression law exponent [-]

% ------------ COMBUSTION DATA ------------
% Oxidizer mass flow, thermochemistry lookup name and O/F - pressure bounds.
oxidizer_lookup = params.combustion.oxidizer_lookup;
mox = params.combustion.mdot_ox;  % Oxidizer mass flow rate [kg/s]
ofmin = params.combustion.of_min; % Lower O/F bound [-]
ofmax = params.combustion.of_max; % Upper O/F bound [-]
pmax = params.combustion.p_max;   % Upper chamber pressure bound [Pa]
pmin = params.combustion.p_min;   % Lower chamber pressure bound [Pa]

% Load the thermochemistry interpolants (temperature, gas constant, gamma).
interp_data = load(fullfile(repoRoot, "LookupTable", oxidizer_lookup, ...
    "interpolant_" + oxidizer_lookup + ".mat"));
T_interp = interp_data.interpolant.T_interp;             % Chamber temperature [K]
R_interp = interp_data.interpolant.R_interp;             % Gas constant [J/kgK]
gamma_interp = interp_data.interpolant.gamma_interp;     % Heat capacity ratio [-]

% Thermochemistry as functions of (O/F, chamber pressure in bar).
T_fun_of_p = @(of,p_bar) T_interp(p_bar, of);
R_fun_of_p = @(of,p_bar) R_interp(p_bar, of);
k_fun_of_p = @(of,p_bar) gamma_interp(p_bar, of);

% Central-difference derivatives of R*T for the fine ODE model.
dp_bar = 1e-3;                    % Pressure step for differentiation [bar]
dof = 1e-3;                       % O/F step for differentiation [-]
RT_fun = @(of,p_bar) R_interp(p_bar, of) * T_interp(p_bar, of);
dRTdp_fun_of_p = @(of,p_bar) (RT_fun(of, p_bar + dp_bar) - RT_fun(of, p_bar - dp_bar)) / (2*dp_bar);
dRTdOF_fun_of_p = @(of,p_bar) (RT_fun(of + dof, p_bar) - RT_fun(of - dof, p_bar)) / (2*dof);

% Time grid, ambient pressure and output / model flags.
tmax = params.time.tmax;          % Simulation end time [s]
time_output = 0:params.time.dt_output:tmax;
pamb = params.engine.pamb;        % Ambient pressure [Pa]
fine_ode_boolean = params.time.fine_ode;
do_plots = params.output.do_plots;
do_animation = params.output.do_animation;
animation_max_frames = params.output.animation_max_frames;
animation_pause = params.output.animation_pause;
animation_autoplay = params.output.animation_autoplay;

% ------------ ENGINE DATA ------------
% Casing, chamber and nozzle dimensions.
ext_diameter = params.engine.ext_diameter;               % Casing external diameter [m]
chamber_length = params.engine.chamber_length;           % Grain / chamber length [m]
throat_diameter = params.engine.throat_diameter;         % Throat diameter [m]
eps = params.engine.eps;                                 % Nozzle area ratio [-]
At = 0.25*pi*(throat_diameter^2);                        % Throat area [m2]

%% STEP 1: build MDF geometry lookup
% Configure and run the Marker-Distance-Field solver, then turn its output
% into a strictly increasing burnback lookup of port area and perimeter.

mdf_cfg = struct();
mdf_cfg.mode = "lookup";
mdf_cfg.geometry_type = geometry_type;
mdf_cfg.meshdata = meshdata;
mdf_cfg.npoints = npoints;
mdf_cfg.casing_radius = ext_diameter / 2;                % [m]
mdf_cfg.t_vec = time_output;
mdf_cfg.h = mdf_cfg.casing_radius / params.mdf.grid_divisions;   % Grid spacing [m]
mdf_cfg.use_bwdist = params.mdf.use_bwdist;
mdf_cfg.use_parallel = params.mdf.use_parallel;
mdf_cfg.do_plots = false;
mdf_cfg.v_reg = 0;
mdf_cfg.geom_opts.use_lookup = true;
mdf_cfg.geom_opts.n_lookup = params.mdf.n_lookup;
mdf_cfg.geom_opts.perimeter_from_area = params.mdf.perimeter_from_area;
mdf_cfg.geom_opts.store_contours = do_animation;
mdf_cfg.geom_opts.smooth_prefix = params.mdf.smooth_prefix_frac;

mdf = burnback_mdf_main(mdf_cfg);
lookup = prepare_geometry_lookup(mdf);

% Initial port area and perimeter at zero burnback.
Ap0 = interp_geometry(0, lookup.b, lookup.Ap);
perim0 = interp_geometry(0, lookup.b, lookup.perim);


%% STEP 2: initialize chamber, st.st.
% Pack the chamber model variables and solve for the steady-state pressure
% that zeroes the chamber balance Z, guarding against a missing sign change.

vars = [];
vars.geometry.port_area = Ap0;
vars.geometry.burning_perimeter = perim0;
vars.geometry.grain_length = chamber_length;
vars.geometry.throat_area = At;
vars.fuel.a_rf = a_rf;
vars.fuel.n_rf = n_rf;
vars.fuel.rho_f = rho_f;
vars.combustion.mdot_ox = mox;
vars.combustion.Tc_fun = T_fun_of_p;
vars.combustion.R_fun = R_fun_of_p;
vars.combustion.k_fun = k_fun_of_p;
vars.combustion.of_min = ofmin;
vars.combustion.of_max = ofmax;
vars.combustion.p_min = pmin;
vars.combustion.p_max = pmax;

% Guard against missing sign change before calling fzero.
[zmin, props_min] = Z_chamber_stst(pmin, vars);
[zmax, props_max] = Z_chamber_stst(pmax, vars);
if zmin == 0
    p0 = pmin;
elseif zmax == 0
    p0 = pmax;
elseif sign(zmin) == sign(zmax)
    if zmin > 0
        error(['No steady-state root in [pmin, pmax]. ' ...
            'Z(pmin) and Z(pmax) are positive, so the steady-state pressure ' ...
            'is above pmax. Try increasing p_max or revisiting mass flow/geometry.']);
    else
        error(['No steady-state root in [pmin, pmax]. ' ...
            'Z(pmin) and Z(pmax) are negative, so the steady-state pressure ' ...
            'is below pmin. Try decreasing p_min or revisiting mass flow/geometry.']);
    end
else
    p0 = fzero(@(pc) Z_chamber_stst(pc, vars), 20e5);
end

% Recover the initial operating point at the steady-state pressure.
[~, properties0] = Z_chamber_stst(p0, vars);
of0 = properties0.of;             % Initial mixture ratio [-]
gox0 = properties0.gox;           % Initial oxidizer mass flux [kg/m2s]
fprintf("Initial values:\n");
fprintf("\tpressure = %.1f bar\n", p0*1e-5);
fprintf("\tGox = %.1f kg/m2s\n", gox0);
fprintf("\tO/F = %.2f\n", of0);

%% STEP 3: integrate pressure and burnback
% Integrate the coupled chamber-pressure / burnback ODE with event-based
% termination on geometry, pressure and O/F limits.

% Hand the fine-ODE derivatives to the model when the fine option is active.
if fine_ode_boolean
    vars.combustion.dRTdp_fun_OF_p = dRTdp_fun_of_p;
    vars.combustion.dRTdOF_fun_OF_p = dRTdOF_fun_of_p;
end

% Initial state [chamber pressure; burnback] and solver options.
y0 = [p0; 0];
options = odeset( ...
    "RelTol", 1e-7, ...
    "AbsTol", [1e-3, 1e-9], ...
    "Events", @(t,y) burnback_events(t, y, lookup, vars));

[time, y_in_time] = ode113( ...
    @(t,y) ode_pressure_burnback(t, y, vars, lookup, fine_ode_boolean), ...
    time_output, y0, options);

in_time.pch = y_in_time(:,1);     % Chamber pressure history [Pa]
in_time.burnback = y_in_time(:,2);% Cumulative burnback history [m]

%% STEP 4: use the results to compute everything else
% March over the integrated history and reconstruct the full set of engine
% performance quantities at every output time.

for j = 1:length(time)
    % Geometry at the current burnback.
    b = in_time.burnback(j);
    Ap = interp_geometry(b, lookup.b, lookup.Ap);
    perim = interp_geometry(b, lookup.b, lookup.perim);
    vars.geometry.port_area = Ap;
    vars.geometry.burning_perimeter = perim;

    % Operating point from the steady-state chamber balance.
    pch = in_time.pch(j);
    pch_bar = pch*1e-5;

    [~, properties] = Z_chamber_stst(pch, vars);
    O_F = properties.of;
    GOX = properties.gox;
    rf = properties.rf;

    % Thermochemistry at the (clamped) operating point.
    [O_F_thermo, pch_bar_thermo] = clamp_thermo_inputs(O_F, pch_bar, vars);
    Tch = T_fun_of_p(O_F_thermo, pch_bar_thermo);
    k = k_fun_of_p(O_F_thermo, pch_bar_thermo);
    R = R_fun_of_p(O_F_thermo, pch_bar_thermo);

    % Nozzle expansion to the exit section.
    Me = supersonic_mach_from_area_ratio(eps, k);
    Te = Tch/(1 + 0.5*(k - 1)*(Me^2));
    pe = pch/((1 + 0.5*(k - 1)*(Me^2))^(k/(k - 1)));
    ae = sqrt(k*R*Te);
    ue = Me*ae;

    % Mass balance, characteristic velocity and nozzle mass flow.
    mfuel = mox/O_F;
    mdot_in = mfuel + mox;
    K2 = k*((2/(k + 1))^((k + 1)/(k - 1)));
    cstar = sqrt(R*Tch/K2);
    mdot_nozzle = pch*At/cstar;

    % Thrust and specific impulse.
    thrust = mdot_nozzle*ue + (pe - pamb)*eps*At;
    I = thrust/mdot_nozzle;
    I_input = thrust/mdot_in;
    Isp = thrust/(mdot_nozzle*9.80665);
    Isp_input = thrust/(mdot_in*9.80665);

    % Store the quantities of interest in time.
    in_time.perim(j) = perim;
    in_time.Ap(j) = Ap;
    in_time.O_F(j) = O_F;
    in_time.GOX(j) = GOX;
    in_time.rf(j) = rf;
    in_time.Tch(j) = Tch;
    in_time.mdot_prop(j) = mdot_in;
    in_time.mdot_in(j) = mdot_in;
    in_time.mdot_nozzle(j) = mdot_nozzle;
    in_time.Me(j) = Me;
    in_time.pe(j) = pe;
    in_time.thrust(j) = thrust;
    in_time.I(j) = I;
    in_time.I_input(j) = I_input;
    in_time.Isp(j) = Isp;
    in_time.Isp_input(j) = Isp_input;
end
toc

I_tot = trapz(time, in_time.thrust);
mdot_tot = in_time.thrust ./ (in_time.Isp .* 9.80665);
m_tot = trapz(time, mdot_tot);
I_avg = I_tot / m_tot;
Isp_avg = I_tot / (m_tot * 9.80665);
fprintf('I_tot: %.2f N*s\n', I_tot);
fprintf('I_avg: %.2f m/s\n', I_avg);
fprintf('Isp_avg: %.2f s\n', Isp_avg);
thrust_avg = I_tot/300;
fprintf('T_avg: %.2f N\n', thrust_avg);


%% plotting
% Performance dashboard, geometry history, initial signed-distance field,
% lookup tables and optional burnback animation.

if do_plots
    % -------- ENGINE PERFORMANCE --------
    figure('Name', 'Hybrid Rocket Engine Performance', 'Color', 'w');

    subplot(2,3,1)
    plot(time, in_time.pch.*1e-5, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Chamber Pressure, $p_{ch}$, bar', 'Interpreter', 'latex')
    title('Chamber Pressure', 'Interpreter', 'latex');

    subplot(2,3,2)
    plot(time, in_time.Tch, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Chamber Temperature, $T_{ch}$, K', 'Interpreter', 'latex')
    title('Chamber Temperature', 'Interpreter', 'latex');

    subplot(2,3,3)
    plot(time, in_time.GOX, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Oxidizer Mass Flux, $G_{ox}$, kg/m$^2$s', 'Interpreter', 'latex')
    title('Oxidizer Mass Flux', 'Interpreter', 'latex');

    subplot(2,3,4)
    plot(time, in_time.O_F, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('O/F', 'Interpreter', 'latex')
    title('Mixture Ratio', 'Interpreter', 'latex');

    subplot(2,3,5)
    plot(time, in_time.thrust, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Thrust, N', 'Interpreter', 'latex')
    title('Thrust', 'Interpreter', 'latex');

    subplot(2,3,6)
    plot(time, in_time.Isp, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Specific Impulse, $I_{sp}$, s', 'Interpreter', 'latex')
    title('Specific Impulse', 'Interpreter', 'latex');

    % -------- CHAMBER PRESSURE --------
    figure('Name', 'Chamber Pressure', 'Color', 'w');
    plot(time, in_time.pch.*1e-5, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Chamber Pressure, $p_{ch}$, bar', 'Interpreter', 'latex')
    title('Chamber Pressure', 'Interpreter', 'latex');

    % ------- CHAMBER TEMPERATURE --------
    figure('Name', 'Chamber Temperature', 'Color', 'w');
    plot(time, in_time.Tch, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Chamber Temperature, $T_{ch}$, K', 'Interpreter', 'latex')
    title('Chamber Temperature', 'Interpreter', 'latex');

    % ------- GOX --------
    figure('Name', 'Oxidizer Mass Flux', 'Color', 'w');
    plot(time, in_time.GOX, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Oxidizer Mass Flux, $G_{ox}$, kg/m$^2$s', 'Interpreter', 'latex')
    title('Oxidizer Mass Flux', 'Interpreter', 'latex');

    % ------ O/F --------
    figure('Name', 'Mixture Ratio', 'Color', 'w');
    plot(time, in_time.O_F, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('O/F', 'Interpreter', 'latex')
    title('Mixture Ratio', 'Interpreter', 'latex');

    % -------- Thrust --------
    figure('Name', 'Thrust', 'Color', 'w');
    plot(time, in_time.thrust, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Thrust, N', 'Interpreter', 'latex')
    title('Thrust', 'Interpreter', 'latex');

    % ------- Specific Impulse --------
    figure('Name', 'Specific Impulse', 'Color', 'w');
    plot(time, in_time.Isp, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Specific Impulse, $I_{sp}$, s', 'Interpreter', 'latex')
    title('Specific Impulse', 'Interpreter', 'latex');

    % -------- IMPULSE --------
    figure('Name', 'Impulse', 'Color', 'w');
    plot(time, in_time.I, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Impulse, I, m/s', 'Interpreter', 'latex')
    title('Impulse', 'Interpreter', 'latex');

    % -------- RF --------
    figure('Name', 'Regression Rate', 'Color', 'w');
    plot(time, in_time.rf, 'LineWidth', 1.5);
    grid on;
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Regression Rate, $r_f$, m/s', 'Interpreter', 'latex')
    title('Regression Rate', 'Interpreter', 'latex');

    % -------- MDF GEOMETRY HISTORY --------
    figure('Name', 'MDF Geometry History', 'Color', 'w');
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    nexttile
    plot(time, in_time.burnback, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Burnback, b, m', 'Interpreter', 'latex')
    title('Cumulative Burnback', 'Interpreter', 'latex');

    nexttile
    plot(time, in_time.Ap, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Port Area, $A_p$, $m^2$', 'Interpreter', 'latex')
    title('Port Area', 'Interpreter', 'latex');

    nexttile
    plot(time, in_time.perim, 'LineWidth', 1.5);
    grid on
    xlabel('Time, t, s', 'Interpreter', 'latex')
    ylabel('Burning Perimeter, $P_b$, m', 'Interpreter', 'latex')
    title('Burning Perimeter', 'Interpreter', 'latex');

    % -------- INITIAL SIGNED DISTANCE FIELD --------
    figure('Name', 'MDF phi0', 'Color', 'w');
    plot_phi0(mdf, ext_diameter/2, throat_diameter/2);

    % -------- GEOMETRY LOOKUP TABLES --------
    figure('Name', 'MDF LUT', 'Color', 'w');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile
    plot(lookup.b, lookup.Ap, 'LineWidth', 1.5);
    grid on
    xlabel('Burnback, $d_b$, m', 'Interpreter', 'latex')
    ylabel('Port Area, $A_p$, $m^2$', 'Interpreter', 'latex')
    title('Port Area LUT', 'Interpreter', 'latex');

    nexttile
    plot(lookup.b, lookup.perim, 'LineWidth', 1.5);
    grid on
    xlabel('Burnback, $d_b$, m', 'Interpreter', 'latex')
    ylabel('Burning Perimeter, $P_b$, m', 'Interpreter', 'latex')
    title('Burning Perimeter LUT', 'Interpreter', 'latex');

    % -------- BURNBACK ANIMATION --------
    if do_animation
        animate_combustion_burnback( ...
            mdf, time, in_time, ext_diameter/2, throat_diameter/2, ...
            mdf_cfg.h, animation_max_frames, animation_pause, animation_autoplay);
    end
end

%% FUNCTIONS

function lookup = prepare_geometry_lookup(mdf)
% prepare_geometry_lookup
% Build a clean, strictly increasing burnback lookup of port area and
% burning perimeter from the raw MDF solver output.
% INPUT
%   mdf    / struct / 1x1   MDF solver output with b/area/perimeter lookups
% OUTPUT
%   lookup / struct / 1x1   Filtered burnback lookup with derivatives

% Pull the raw burnback, port area and perimeter columns.
lookup.b = mdf.b_lookup(:);
lookup.Ap = mdf.port_area_lookup(:);
lookup.perim = mdf.perimeter_lookup(:);

% Keep only finite, physically meaningful samples.
valid = isfinite(lookup.b) & isfinite(lookup.Ap) & isfinite(lookup.perim) & ...
    lookup.Ap > 0 & lookup.perim >= 0;
lookup.b = lookup.b(valid);
lookup.Ap = lookup.Ap(valid);
lookup.perim = lookup.perim(valid);

% Enforce uniqueness and a strictly increasing burnback axis.
[lookup.b, unique_idx] = unique(lookup.b, 'stable');
lookup.Ap = lookup.Ap(unique_idx);
lookup.perim = lookup.perim(unique_idx);

if numel(lookup.b) < 2
    error('MDF lookup must contain at least two valid burnback levels.');
end
if any(diff(lookup.b) <= 0)
    error('MDF burnback lookup must be strictly increasing.');
end

% Rebuild port area by integrating the perimeter (dAp/db = perimeter)
% and clamp it inside the raw bounds; store the geometry derivatives.
lookup.Ap_raw = lookup.Ap;
lookup.Ap = lookup.Ap(1) + cumtrapz(lookup.b, lookup.perim);
lookup.Ap = min(max(lookup.Ap, lookup.Ap_raw(1)), max(lookup.Ap_raw));
lookup.dAp_db = lookup.perim;
lookup.dperim_db = gradient(lookup.perim, lookup.b);
lookup.b_max = lookup.b(end);
end

function dy = ode_pressure_burnback(~, y, vars, lookup, fine_ode_boolean)
% ode_pressure_burnback
% Right-hand side of the coupled chamber-pressure / burnback ODE.
% INPUT
%   ~                / double / 1x1   Time (unused)
%   y                / double / 2x1   State [chamber pressure; burnback]
%   vars             / struct / 1x1   Chamber model variables
%   lookup           / struct / 1x1   Geometry burnback lookup
%   fine_ode_boolean / logical/ 1x1   Use the fine pressure model when true
% OUTPUT
%   dy               / double / 2x1   State derivative [dp/dt; regression rate]

pch = y(1);
b = y(2);

% Geometry and its derivatives at the current burnback.
Ap = interp_geometry(b, lookup.b, lookup.Ap);
perim = interp_geometry(b, lookup.b, lookup.perim);
dAp_db = interp_geometry(b, lookup.b, lookup.dAp_db);
dperim_db = interp_geometry(b, lookup.b, lookup.dperim_db);

% Unpack the model variables.
L = vars.geometry.grain_length;
At = vars.geometry.throat_area;
a_rf = vars.fuel.a_rf;
n_rf = vars.fuel.n_rf;
rho_f = vars.fuel.rho_f;
mdot_ox = vars.combustion.mdot_ox;
Tc_fun = vars.combustion.Tc_fun;
R_fun = vars.combustion.R_fun;
k_fun = vars.combustion.k_fun;

% Oxidizer flux, regression rate and mass balance.
Gox = mdot_ox / Ap;
rf = a_rf * (Gox^n_rf);
Ab = perim * L;
mdot_f = rho_f * Ab * rf;
mdot_in = mdot_f + mdot_ox;
O_F = mdot_ox / mdot_f;

% Thermochemistry and nozzle outflow at the (clamped) operating point.
pc_bar = pch * 1e-5;
[O_F_thermo, pc_bar_thermo] = clamp_thermo_inputs(O_F, pc_bar, vars);
Tc = Tc_fun(O_F_thermo, pc_bar_thermo);
R = R_fun(O_F_thermo, pc_bar_thermo);
k = k_fun(O_F_thermo, pc_bar_thermo);
K2 = k*((2/(k + 1))^((k + 1)/(k - 1)));
cstar = sqrt(R*Tc/K2);
mdot_out = pch*At/cstar;
rho_g = pch/(R*Tc);

% Geometry terms fed to the chamber-pressure model.
geometry.grain_length = L;
geometry.port_area = Ap;
geometry.d_area_area = (dAp_db / Ap) * rf;

% Fine model adds the R*T sensitivity and perimeter-growth terms.
if fine_ode_boolean
    dRTdp = vars.combustion.dRTdp_fun_OF_p(O_F_thermo, pc_bar_thermo);
    dRTdOF = vars.combustion.dRTdOF_fun_OF_p(O_F_thermo, pc_bar_thermo);

    fine_ode.deltapc = pc_bar_thermo * dRTdp / (R*Tc);
    fine_ode.deltaOF = O_F_thermo * dRTdOF / (R*Tc);
    fine_ode.n_rf = n_rf;
    fine_ode.d_perim_perim = (dperim_db / max(perim, realmin)) * rf;

    dp_p = ode_chamber_pressure(mdot_in, mdot_out, rho_g, geometry, fine_ode);
else
    dp_p = ode_chamber_pressure(mdot_in, mdot_out, rho_g, geometry);
end

% dp/dt = (relative pressure rate)*pch; db/dt = regression rate.
dy = [dp_p*pch; rf];
end

function [value, isterminal, direction] = burnback_events(~, y, lookup, vars)
% burnback_events
% Event function stopping the integration at geometry, pressure or O/F limits.
% INPUT
%   ~          / double / 1x1   Time (unused)
%   y          / double / 2x1   State [chamber pressure; burnback]
%   lookup     / struct / 1x1   Geometry burnback lookup
%   vars       / struct / 1x1   Chamber model variables
% OUTPUT
%   value      / double / 6x1   Event functions (zero crossing triggers stop)
%   isterminal / double / 6x1   Terminal flags
%   direction  / double / 6x1   Crossing direction per event

pch = y(1);
b = y(2);
perim = interp_geometry(b, lookup.b, lookup.perim);
O_F = mixture_ratio_from_geometry(b, lookup, vars);

value = [
    lookup.b_max - b
    perim
    pch - vars.combustion.p_min
    vars.combustion.p_max - pch
    O_F - vars.combustion.of_min
    vars.combustion.of_max - O_F
];
isterminal = ones(size(value));
direction = [-1; -1; -1; -1; 0; 0];
end

function O_F = mixture_ratio_from_geometry(b, lookup, vars)
% mixture_ratio_from_geometry
% Mixture ratio O/F implied by the grain geometry at a given burnback.
% INPUT
%   b      / double / 1x1   Burnback [m]
%   lookup / struct / 1x1   Geometry burnback lookup
%   vars   / struct / 1x1   Chamber model variables
% OUTPUT
%   O_F    / double / 1x1   Mixture ratio [-]

Ap = interp_geometry(b, lookup.b, lookup.Ap);
perim = interp_geometry(b, lookup.b, lookup.perim);
Gox = vars.combustion.mdot_ox / Ap;
rf = vars.fuel.a_rf * (Gox^vars.fuel.n_rf);
mdot_f = vars.fuel.rho_f * perim * vars.geometry.grain_length * rf;
O_F = vars.combustion.mdot_ox / mdot_f;
end

function [O_F_eval, p_bar_eval] = clamp_thermo_inputs(O_F, p_bar, vars)
% clamp_thermo_inputs
% Clamp O/F and pressure to the lookup bounds before sampling thermochemistry.
% INPUT
%   O_F        / double / 1x1   Mixture ratio [-]
%   p_bar      / double / 1x1   Chamber pressure [bar]
%   vars       / struct / 1x1   Chamber model variables (holds the bounds)
% OUTPUT
%   O_F_eval   / double / 1x1   Clamped mixture ratio [-]
%   p_bar_eval / double / 1x1   Clamped chamber pressure [bar]

O_F_eval = min(max(O_F, vars.combustion.of_min), vars.combustion.of_max);
p_bar_eval = min(max(p_bar, vars.combustion.p_min*1e-5), vars.combustion.p_max*1e-5);
end

function value = interp_geometry(b, b_lookup, values)
% interp_geometry
% Linear interpolation of a geometry quantity at a clamped burnback.
% INPUT
%   b        / double / 1x1   Query burnback [m]
%   b_lookup / double / Nx1   Burnback axis [m]
%   values   / double / Nx1   Quantity sampled on the burnback axis
% OUTPUT
%   value    / double / 1x1   Interpolated quantity

b = min(max(b, b_lookup(1)), b_lookup(end));
value = interp1(b_lookup, values, b, "linear");
end

function Me = supersonic_mach_from_area_ratio(area_ratio, k)
% supersonic_mach_from_area_ratio
% Supersonic exit Mach number for a given nozzle area ratio and gamma.
% INPUT
%   area_ratio / double / 1x1   Nozzle area ratio Ae/At [-]
%   k          / double / 1x1   Heat capacity ratio [-]
% OUTPUT
%   Me         / double / 1x1   Supersonic exit Mach number [-]

if area_ratio < 1
    error('Nozzle area ratio must be >= 1.');
end

% Bracket the supersonic branch, then solve the area-Mach relation.
lo = 1 + 1e-8;
hi = 2;
while area_mach_ratio(hi, k) < area_ratio
    hi = hi * 1.5;
    if hi > 100
        error('Could not bracket supersonic exit Mach number.');
    end
end

Me = fzero(@(M) area_mach_ratio(M, k) - area_ratio, [lo, hi]);
end

function eps_val = area_mach_ratio(M, k)
% area_mach_ratio
% Isentropic area ratio A/A* as a function of Mach number and gamma.
% INPUT
%   M       / double / 1x1   Mach number [-]
%   k       / double / 1x1   Heat capacity ratio [-]
% OUTPUT
%   eps_val / double / 1x1   Area ratio A/A* [-]

eps_val = 1/M * sqrt(((1 + 0.5*(k - 1)*(M^2))/(0.5*(k + 1)))^((k + 1)/(k - 1)));
end

function plot_phi0(mdf, casing_radius, throat_radius)
% plot_phi0
% Plot the initial signed distance field phi0 with case and throat outlines.
% INPUT
%   mdf           / struct / 1x1   MDF solver output (phi0 and grid)
%   casing_radius / double / 1x1   Case radius [m] (optional)
%   throat_radius / double / 1x1   Throat radius [m] (optional)
% OUTPUT
%   None

phi0 = mdf.phi0;
X = mdf.X;
Y = mdf.Y;
xv = mdf.xvec;
yv = mdf.yvec;

% Rebuild the grid if the solver did not store the meshgrid.
if isempty(X) || isempty(Y)
    [X, Y] = meshgrid(xv, yv);
end

% Filled field with the zero level set highlighted.
contourf(X, Y, phi0, 50, 'LineColor', 'none');
hold on;
[~, contour_handle] = contour(X, Y, phi0, [0 0], 'k-', 'LineWidth', 2.0);
set(contour_handle, 'DisplayName', 'phi0 = 0');
legend_handles = contour_handle;

% Overlay case and throat half circles when provided.
if nargin >= 2 && ~isempty(casing_radius)
    case_handle = plot_half_circle(casing_radius, 'k-', 'case');
    legend_handles(end+1) = case_handle;
end
if nargin >= 3 && ~isempty(throat_radius)
    throat_handle = plot_half_circle(throat_radius, 'k--', 'throat');
    legend_handles(end+1) = throat_handle;
end
cb = colorbar;
cb.Label.String = 'phi0 [m]';
colormap(gca, 'turbo');
xlabel('x, m', 'Interpreter', 'latex');
ylabel('y, m', 'Interpreter', 'latex');
title('phi0(x,y) - Initial Signed Distance Field (half domain)', 'Interpreter', 'latex');
axis equal;
xlim([min(xv), max(xv)]);
ylim([min(yv), max(yv)]);
grid on;
legend(legend_handles, 'Location', 'northwest');
end

function handle = plot_half_circle(radius, style, display_name)
% plot_half_circle
% Plot the upper half of a circle of given radius.
% INPUT
%   radius       / double / 1x1   Circle radius [m]
%   style        / char   / 1xN   Line style spec
%   display_name / char   / 1xN   Legend entry
% OUTPUT
%   handle       / Line   / 1x1   Handle to the plotted half circle

theta = linspace(0, pi, 500);
handle = plot(radius*cos(theta), radius*sin(theta), style, 'DisplayName', display_name);
end

function animate_combustion_burnback(mdf, time, in_time, casing_radius, throat_radius, h, max_frames, pause_time, autoplay)
% animate_combustion_burnback
% Interactive animation of the burning surface together with the perimeter
% and port-area history (play button + slider).
% INPUT
%   mdf           / struct / 1x1   MDF solver output (contours lookup)
%   time          / double / Nx1   Output time vector [s]
%   in_time       / struct / 1x1   Time histories (burnback, perim, Ap)
%   casing_radius / double / 1x1   Case radius [m]
%   throat_radius / double / 1x1   Throat radius [m]
%   h             / double / 1x1   Grid spacing [m]
%   max_frames    / double / 1x1   Max number of frames (optional)
%   pause_time    / double / 1x1   Pause between frames [s] (optional)
%   autoplay      / logical/ 1x1   Auto-start playback (optional)
% OUTPUT
%   None

% Default the optional arguments.
if nargin < 7 || isempty(max_frames)
    max_frames = 120;
end
if nargin < 8 || isempty(pause_time)
    pause_time = 0.03;
end
if nargin < 9 || isempty(autoplay)
    autoplay = false;
end

% Select the frame indices to display and the burning cutoff radius.
time = time(:);
n_time = numel(time);
if n_time == 0
    return
end

frame_idx = unique(round(linspace(1, n_time, min(max_frames, n_time))), 'stable');
n_frames = numel(frame_idx);
burn_limit = max(casing_radius - 0.5*h, 0);

% -------- FIGURE AND GEOMETRY AXES --------
fig = figure('Name', 'Combustion Burnback Animation', 'Color', 'w');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

geometry_axes = nexttile(layout, 1);
hold(geometry_axes, 'on');
axis(geometry_axes, 'equal');
grid(geometry_axes, 'on');
box(geometry_axes, 'on');
xlabel(geometry_axes, 'x, m');
ylabel(geometry_axes, 'y, m');
xlim(geometry_axes, casing_radius * [-1.08, 1.08]);
ylim(geometry_axes, casing_radius * [-1.08, 1.08]);
plot_circle_on_axes(geometry_axes, casing_radius, 'k-', 'case');
plot_circle_on_axes(geometry_axes, burn_limit, 'k:', 'burning cutoff');
plot_circle_on_axes(geometry_axes, throat_radius, 'k--', 'throat');

% Animated burning surface and case-contact (non-burning) line.
burning_line = plot(geometry_axes, nan, nan, '-', ...
    'Color', [0.00 0.45 0.74], 'LineWidth', 1.7, 'DisplayName', 'burning surface');
wall_line = plot(geometry_axes, nan, nan, '-', ...
    'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, 'DisplayName', 'case contact, no burn');
title_handle = title(geometry_axes, '');
legend(geometry_axes, 'Location', 'best');

% -------- HISTORY AXES (PERIMETER AND PORT AREA) --------
history_axes = nexttile(layout, 2);
hold(history_axes, 'on');
grid(history_axes, 'on');
box(history_axes, 'on');
xlabel(history_axes, 'time, t, s');
ylabel(history_axes, 'P_b, m');
plot(history_axes, time, in_time.perim, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, ...
    'DisplayName', 'P_b');
perim_marker = plot(history_axes, time(1), in_time.perim(1), 'o', ...
    'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], ...
    'DisplayName', 'current');
yyaxis(history_axes, 'right');
ylabel(history_axes, 'A_p, m^2');
plot(history_axes, time, in_time.Ap, 'Color', [0.00 0.45 0.74], 'LineWidth', 1.4, ...
    'DisplayName', 'A_p');
area_marker = plot(history_axes, time(1), in_time.Ap(1), 's', ...
    'Color', [0.00 0.45 0.74], 'MarkerFaceColor', [0.00 0.45 0.74], ...
    'DisplayName', 'current');
yyaxis(history_axes, 'left');

% -------- PLAYBACK CONTROLS --------
if n_frames > 1
    slider_step = [1/(n_frames - 1), min(1, 10/(n_frames - 1))];
else
    slider_step = [1, 1];
end

uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.44 0.02 0.10 0.04], 'String', 'Play', ...
    'Callback', @play_animation);
frame_slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.20 0.02 0.22 0.04], 'Min', 1, 'Max', n_frames, ...
    'Value', 1, 'SliderStep', slider_step, 'Callback', @slider_changed);
frame_text = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.56 0.02 0.24 0.04], 'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left');

% Render the first frame and optionally start playback.
update_frame(1);
if autoplay
    play_animation([], []);
end

    function slider_changed(source, ~)
        % Jump to the frame selected on the slider.
        update_frame(round(get(source, 'Value')));
    end

    function play_animation(~, ~)
        % Play from the current slider position to the last frame.
        start_frame = round(get(frame_slider, 'Value'));
        for frame_number = start_frame:n_frames
            if ~ishandle(fig)
                return
            end
            update_frame(frame_number);
            pause(pause_time);
        end
    end

    function update_frame(frame_number)
        % Redraw the burning surface and history markers for a given frame.
        frame_number = max(1, min(n_frames, round(frame_number)));
        k = frame_idx(frame_number);

        b = in_time.burnback(k);
        polys = contours_at_burnback(mdf, b);
        [xb, yb, xw, yw] = split_burning_and_wall_segments(polys, casing_radius, h);

        set(burning_line, 'XData', xb, 'YData', yb);
        set(wall_line, 'XData', xw, 'YData', yw);
        set(perim_marker, 'XData', time(k), 'YData', in_time.perim(k));
        yyaxis(history_axes, 'right');
        set(area_marker, 'XData', time(k), 'YData', in_time.Ap(k));
        yyaxis(history_axes, 'left');
        set(frame_slider, 'Value', frame_number);

        set(title_handle, 'String', sprintf( ...
            't = %.2f s, b = %.4f m, P_b = %.4f m', ...
            time(k), b, in_time.perim(k)));
        set(frame_text, 'String', sprintf('Frame %d/%d', frame_number, n_frames));
        drawnow limitrate;
    end
end

function polys = contours_at_burnback(mdf, b)
% contours_at_burnback
% Return the stored half-contour set closest to a given burnback level.
% INPUT
%   mdf   / struct / 1x1   MDF solver output (contours lookup)
%   b     / double / 1x1   Query burnback [m]
% OUTPUT
%   polys / cell   / 1xP   Half-domain contour polylines

[~, idx] = min(abs(mdf.b_lookup(:) - b));
polys = mdf.half_contours_lookup{idx};
end

function [xb, yb, xw, yw] = split_burning_and_wall_segments(polys, casing_radius, h)
% split_burning_and_wall_segments
% Split a half-contour set (mirrored about x) into burning and case-contact
% segments, NaN-separated for a single plot call.
% INPUT
%   polys         / cell   / 1xP   Half-domain contour polylines
%   casing_radius / double / 1x1   Case radius [m]
%   h             / double / 1x1   Grid spacing [m]
% OUTPUT
%   xb, yb        / double / Mx1   Burning segment coordinates [m]
%   xw, yw        / double / Mx1   Case-contact segment coordinates [m]

burn_limit2 = max(casing_radius - 0.5*h, 0)^2;
xb = [];
yb = [];
xw = [];
yw = [];

for p = 1:numel(polys)
    xy_top = polys{p};
    if isempty(xy_top) || size(xy_top, 1) < 2
        continue
    end

    % Top half of the contour.
    [x_b, y_b, x_w, y_w] = split_polyline_segments(xy_top, burn_limit2);

    xb = [xb; x_b; nan]; %#ok<AGROW>
    yb = [yb; y_b; nan]; %#ok<AGROW>
    xw = [xw; x_w; nan]; %#ok<AGROW>
    yw = [yw; y_w; nan]; %#ok<AGROW>

    % Mirrored bottom half of the contour.
    xy_bot = [xy_top(:,1), -xy_top(:,2)];
    [x_b, y_b, x_w, y_w] = split_polyline_segments(xy_bot, burn_limit2);

    xb = [xb; x_b; nan]; %#ok<AGROW>
    yb = [yb; y_b; nan]; %#ok<AGROW>
    xw = [xw; x_w; nan]; %#ok<AGROW>
    yw = [yw; y_w; nan]; %#ok<AGROW>
end
end

function [x_b, y_b, x_w, y_w] = split_polyline_segments(xy, burn_limit2)
% split_polyline_segments
% Classify each polyline edge as burning or case-contact by its midpoint radius.
% INPUT
%   xy          / double / Nx2   Polyline coordinates [m]
%   burn_limit2 / double / 1x1   Squared burning cutoff radius [m2]
% OUTPUT
%   x_b, y_b    / double / Mx1   Burning edge coordinates [m]
%   x_w, y_w    / double / Mx1   Case-contact edge coordinates [m]

dx = diff(xy(:,1));
dy = diff(xy(:,2));
xm = xy(1:end-1,1) + 0.5*dx;
ym = xy(1:end-1,2) + 0.5*dy;
burning = (xm.^2 + ym.^2) <= burn_limit2;

[x_b, y_b] = segment_nan_data(xy, burning);
[x_w, y_w] = segment_nan_data(xy, ~burning);
end

function [x_plot, y_plot] = segment_nan_data(xy, mask)
% segment_nan_data
% Expand the masked polyline edges into NaN-separated drawable segments.
% INPUT
%   xy     / double  / Nx2   Polyline coordinates [m]
%   mask   / logical / (N-1)x1   Edge selection mask
% OUTPUT
%   x_plot / double  / Mx1   Segment x coordinates, NaN-separated [m]
%   y_plot / double  / Mx1   Segment y coordinates, NaN-separated [m]

x_plot = nan(3*nnz(mask), 1);
y_plot = nan(3*nnz(mask), 1);
out = 1;
for i = find(mask(:))'
    x_plot(out:out+2) = [xy(i,1); xy(i+1,1); nan];
    y_plot(out:out+2) = [xy(i,2); xy(i+1,2); nan];
    out = out + 3;
end
end

function plot_circle_on_axes(ax, radius, style, display_name)
% plot_circle_on_axes
% Plot a full circle of given radius on a specific axes handle.
% INPUT
%   ax           / Axes   / 1x1   Target axes
%   radius       / double / 1x1   Circle radius [m]
%   style        / char   / 1xN   Line style spec
%   display_name / char   / 1xN   Legend entry
% OUTPUT
%   None

theta = linspace(0, 2*pi, 500);
plot(ax, radius*cos(theta), radius*sin(theta), style, "DisplayName", display_name);
end
