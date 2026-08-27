function [Z, properties] = Z_chamber_stst(pc,vars)

% define error for mass balance in combustion chamber 
% in steady state condition

%% extract data
% geometry
Ap = vars.geometry.port_area;
perimB = vars.geometry.burning_perimeter;
L = vars.geometry.grain_length;
At = vars.geometry.throat_area;
% fuel
a_rf = vars.fuel.a_rf;
n_rf = vars.fuel.n_rf;
rho_f = vars.fuel.rho_f;
% combustion
mdot_ox = vars.combustion.mdot_ox;
Tc_fun = vars.combustion.Tc_fun;
R_fun = vars.combustion.R_fun;
k_fun = vars.combustion.k_fun;

%% calculate
% compute mdotin
Gox = mdot_ox/Ap;
rf = a_rf*(Gox^n_rf);
Ab = perimB*L;
mdot_f = rho_f*Ab*rf;
mdot_in = mdot_f + mdot_ox; % portata che entra nel sistema

% compute mdotout
O_F = mdot_ox/mdot_f;
pc_bar = pc*1e-5;
validate_thermo_inputs(O_F, pc_bar, vars);
Tc = Tc_fun(O_F,pc_bar); % funzioni che variano in base all'ossidante
R = R_fun(O_F,pc_bar);
k = k_fun(O_F,pc_bar);
K2 = k*((2/(k+1))^((k+1)/(k-1))); % squared Vandenkerckhove
cstar = sqrt(R*Tc/K2);
mdot_out = pc*At/cstar;

% error
Z = mdot_in-mdot_out;

% properties
properties.of = O_F;
properties.gox = Gox;
properties.rf = rf;

end

function validate_thermo_inputs(O_F, pc_bar, vars)
if isfield(vars.combustion, 'of_min') && isfield(vars.combustion, 'of_max')
    if O_F < vars.combustion.of_min || O_F > vars.combustion.of_max
        error('O/F %.6g outside thermochemical lookup range [%.6g, %.6g].', ...
            O_F, vars.combustion.of_min, vars.combustion.of_max);
    end
end

if isfield(vars.combustion, 'p_min') && isfield(vars.combustion, 'p_max')
    p_min_bar = vars.combustion.p_min * 1e-5;
    p_max_bar = vars.combustion.p_max * 1e-5;
    if pc_bar < p_min_bar || pc_bar > p_max_bar
        error('Pressure %.6g bar outside thermochemical lookup range [%.6g, %.6g] bar.', ...
            pc_bar, p_min_bar, p_max_bar);
    end
end
end
