clear; clc; close all;

%% Run CEA to build lookup tables

% Rocket problem, infinite combustion chamber area, frozen equilibrium
%   - Oxidizer: O2 (300 K), O2(L) (90 K), H2O2 (300 K), N2O (300 K), N2O4 (300 K)
%   - Fuel: HTPB, C 7.075, H 10.650, O 0.223, N 0.063, DeltaH0_f = -58.0 kJ/mol (300 K)
% Required parameters for the lookup table:
%   - Flame temperature, Tf         [K]
%   - Molar mass, Mm                [g/mol]
%   - Specific heat ratio, gamma    [-]

p_vec = [1, 1.01325, 2:1:100];      % [bar]

lookupTables = repmat(struct( ...
    'oxidizer', '', ...
    'components', [], ...
    'oxidizerTemp', [], ...
    'of_vec', [], ...
    'T', [], ...
    'R', [], ...
    'gamma', []), 1, 8);

lookupTables(1).oxidizer = 'O2';
lookupTables(1).components = {'O2', 100};
lookupTables(1).oxidizerTemp = 300;
lookupTables(1).of_vec = [1:0.1:3, 3.5:0.5:20];
lookupTables(2).oxidizer = 'O2(L)';
lookupTables(2).components = {'O2(L)', 100};
lookupTables(2).oxidizerTemp = 90;
lookupTables(2).of_vec = [1:0.1:3, 3.5:0.5:20];                 % since optimum o/f for o2 circa 2
lookupTables(3).oxidizer = 'H2O2';
lookupTables(3).components = {'H2O2', 100};
lookupTables(3).oxidizerTemp = 300;
lookupTables(3).of_vec = [1:0.5:4, 4.5:0.1:6.5, 7:0.5:20];      % since optimum o/f for pure h2o2 circa 5.5
lookupTables(4).oxidizer = 'H2O2_90';
lookupTables(4).components = {'H2O2', 90; 'H2O', 10};
lookupTables(4).oxidizerTemp = 300;
lookupTables(4).of_vec = [1:0.5:4, 4.5:0.1:6.5, 7:0.5:20];      % similar range to pure H2O2
lookupTables(5).oxidizer = 'N2O';
lookupTables(5).components = {'N2O', 100};
lookupTables(5).oxidizerTemp = 300;
lookupTables(5).of_vec = [1:0.5:5.5, 6:0.1:8, 8.5:0.5:20];      % since optimum o/f for n2o circa 7
lookupTables(6).oxidizer = 'N2O4';
lookupTables(6).components = {'N2O4', 100};
lookupTables(6).oxidizerTemp = 300;
lookupTables(6).of_vec = [1:0.5:1.5, 2:0.1:4.5, 5:0.5:20];      % since optimum o/f for n2o circa 3.5
lookupTables(7).oxidizer = 'IRFNA';
lookupTables(7).components = {'IRFNA', 100};
lookupTables(7).oxidizerTemp = 300;
lookupTables(7).of_vec = [1:0.5:3, 3.5:0.1:6.5, 7:0.5:20];      % since optimum o/f for IRFNA circa 4.3
lookupTables(8).oxidizer = 'FLOX_2F2_1O2';
lookupTables(8).components = {'F2(L)', 70.37; 'O2(L)', 29.63};
lookupTables(8).oxidizerTemp = 87;
lookupTables(8).of_vec = [1:0.5:2, 2.5:0.1:6, 6.5:0.5:20];     % since optimum o/f for FLOX circa 3.3

nOx = numel(lookupTables);

for oxIdx = 1:nOx
    lookupTables(oxIdx) = buildLookupTableForOxidizer(p_vec, lookupTables(oxIdx));
end

%% Export the lookup tables (MAT and CSV, per oxidizer)

% BE CAREFUL IF YOU ARE RUNNING THIS FROM VS CODE, THE SCRIPT DIRECTORY
% MIGHT POINT TO A RANDOM TEMP FOLDER IN APPDATA
scriptDir = fileparts(mfilename('fullpath'));

save(fullfile(scriptDir, 'lookupTable.mat'), 'p_vec', 'lookupTables');

for oxIdx = 1:nOx
    exportLookupTableCsvs(scriptDir, p_vec, lookupTables(oxIdx));
end

%% Visualize the lookup tables (per oxidizer)

for oxIdx = 1:nOx
    visualizeLookupTable(p_vec, lookupTables(oxIdx));
end

%% Build and save interpolants

interpolants = repmat(struct('oxidizer','', 'of_vec', [], 'T_interp', [], 'R_interp', [], 'gamma_interp', []), 1, nOx);

for oxIdx = 1:nOx
    lt = lookupTables(oxIdx);
    interpolants(oxIdx).oxidizer = lt.oxidizer;
    interpolants(oxIdx).of_vec = lt.of_vec;
    interpolants(oxIdx).T_interp = griddedInterpolant({p_vec, lt.of_vec}, lt.T, 'linear', 'none');
    interpolants(oxIdx).R_interp = griddedInterpolant({p_vec, lt.of_vec}, lt.R, 'linear', 'none');
    interpolants(oxIdx).gamma_interp = griddedInterpolant({p_vec, lt.of_vec}, lt.gamma, 'linear', 'none');

    safeName = regexprep(lt.oxidizer, '[^a-zA-Z0-9()]', '_');
    outDir = fullfile(scriptDir, safeName);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    interpolant = interpolants(oxIdx);
    save(fullfile(outDir, sprintf('interpolant_%s.mat', safeName)), 'p_vec', 'interpolant');
end

save(fullfile(scriptDir, 'lookupTableInterpolants.mat'), 'p_vec', 'interpolants');

%% Functions

function lookupTable = buildLookupTableForOxidizer(p_vec, lookupTable)

of_vec = lookupTable.of_vec;
nOf = numel(of_vec);

lookupTable.T = zeros(numel(p_vec), nOf);
lookupTable.R = zeros(numel(p_vec), nOf);
lookupTable.gamma = zeros(numel(p_vec), nOf);

fprintf('\nCalculating for oxidizer: %s (T=%.0f K)\n', lookupTable.oxidizer, lookupTable.oxidizerTemp);
for ii = 1:numel(p_vec)
    for jj = 1:nOf
        try
        p = p_vec(ii);
        of = of_vec(jj);
        %fprintf('Progress: %s, %d/%d (p=%.2f bar, O/F=%.2f)\n', lookupTable.oxidizer, ii, numel(p_vec), p, of);
        
        % Build oxidizer specification from components
        oxid_spec = {};
        for cc = 1:size(lookupTable.components, 1)
            comp_name = lookupTable.components{cc, 1};
            comp_wt = lookupTable.components{cc, 2};
            oxid_spec = [oxid_spec, 'oxid', comp_name, 'wt%', comp_wt, 't,k', lookupTable.oxidizerTemp];
        end
        
        CEA_output = CEA( ...
            'problem', 'o/f', of, 'rocket', 'frozen', 'nfz', 1, 'p,bar', p, ...
            'reactants', ...
                oxid_spec{:}, ...
                'fuel', 'HTPB', 'wt%', 100, 't,k', 300, 'h,kj/mol', -58, 'C', 7.075, 'H', 10.650, 'O', 0.223, 'N', 0.063, ...
            'output', 'massf', 'short', 'transport', ...
            'end');

        Tf = CEA_output.output.froz.temperature(1);      % [K]
        Mm = CEA_output.output.froz.mw(1);               % [g/mol]
        cp = CEA_output.output.froz.cp_tran.froz(1);     % [kJ/kgK]
        cp = cp * 1e3;                                   % [J/kgK]
        R = 8314.5 / Mm;                                 % [J/kgK]
        cv = cp - R;                                     % [J/kgK]
        gamma = cp / cv;                                 % [-]

        lookupTable.T(ii, jj) = Tf;
        lookupTable.R(ii, jj) = R;
        lookupTable.gamma(ii, jj) = gamma;
        catch
        end
    end
end

end

function exportLookupTableCsvs(scriptDir, p_vec, lookupTable)

    safeName = regexprep(lookupTable.oxidizer, '[^a-zA-Z0-9()]', '_');
    variableNames = matlab.lang.makeValidName("of_" + string(lookupTable.of_vec));

    T_table = array2table(lookupTable.T, 'VariableNames', variableNames);
    T_table = addvars(T_table, p_vec(:), 'Before', 1, 'NewVariableNames', 'p_bar');

    R_table = array2table(lookupTable.R, 'VariableNames', variableNames);
    R_table = addvars(R_table, p_vec(:), 'Before', 1, 'NewVariableNames', 'p_bar');

    g_table = array2table(lookupTable.gamma, 'VariableNames', variableNames);
    g_table = addvars(g_table, p_vec(:), 'Before', 1, 'NewVariableNames', 'p_bar');

    outDir = fullfile(scriptDir, safeName);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    writetable(T_table, fullfile(outDir, sprintf('T_lookupTable_%s.csv', safeName)));
    writetable(R_table, fullfile(outDir, sprintf('R_lookupTable_%s.csv', safeName)));
    writetable(g_table, fullfile(outDir, sprintf('gamma_lookupTable_%s.csv', safeName)));

    lookup = lookupTable;
    save(fullfile(outDir, sprintf('lookupTable_%s.mat', safeName)), 'p_vec', 'lookup');

end

function visualizeLookupTable(p_vec, lookupTable)

figure('Name', sprintf('Lookup - %s', lookupTable.oxidizer));
subplot(1,3,1);
surf(lookupTable.of_vec, p_vec, lookupTable.T, 'EdgeColor', 'none');
xlabel('O/F Ratio');
ylabel('Pressure (bar)');
zlabel('Flame Temperature (K)');
title(sprintf('T - %s', lookupTable.oxidizer));

subplot(1,3,2);
surf(lookupTable.of_vec, p_vec, lookupTable.R, 'EdgeColor', 'none');
xlabel('O/F Ratio');
ylabel('Pressure (bar)');
zlabel('Gas Constant (J/kgK)');
title(sprintf('R - %s', lookupTable.oxidizer));

subplot(1,3,3);
surf(lookupTable.of_vec, p_vec, lookupTable.gamma, 'EdgeColor', 'none');
xlabel('O/F Ratio');
ylabel('Pressure (bar)');
zlabel('Specific Heat Ratio (-)');
title(sprintf('gamma - %s', lookupTable.oxidizer));

end