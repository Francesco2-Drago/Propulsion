function dp_p = ode_chamber_pressure(mdot_in,mdot_out,rho_g,geometry,fine_ode)

% derivazione matematica su pdf per ottenere queste formule
% compute normalized pressure variation in time in chamber

%% get data
L = geometry.grain_length;
Ap = geometry.port_area;
d_area_area = geometry.d_area_area;
% due casi per considerare o meno i delta finali della derivazione matemat
if nargin < 5
    n_rf = 0;
    d_perim_perim = 0;
    deltapc = 0;
    deltaOF = 0;
else
    n_rf = fine_ode.n_rf;
    d_perim_perim = fine_ode.d_perim_perim;
    deltapc = fine_ode.deltapc;
    deltaOF = fine_ode.deltaOF;
end

term_mass_balance = (mdot_in-mdot_out) / (L*Ap*rho_g);
term_of_perimeter = deltaOF * d_perim_perim;
term_area_coupling = (1 - n_rf * deltaOF) * d_area_area;

dp_p = (1 / (1 - deltapc)) * (term_mass_balance - term_of_perimeter - term_area_coupling);

end

