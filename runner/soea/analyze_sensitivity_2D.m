% Parallel
p = gcp('nocreate');
if isempty(p)
    parpool;
end

% Common control
clear;
close all;

% Experiment settings
solver = 'SPS_L_SHADE_EIG';
fitfun = 'cec15_f1';
D = 30;
maxfunevals = D * 1e4;
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);

% Parameter baseline
solverOptions.NP = 75 + 345;
solverOptions.F = 0.4076;
solverOptions.CR = 0.6209;
solverOptions.ER = 0.1399;
solverOptions.p = 0.1340;
solverOptions.H = 90;
solverOptions.Q = 194;
solverOptions.Ar = 1.9630;
solverOptions.cw = 0.0581;
solverOptions.erw = 0.6807;
solverOptions.CRmin = 0.3046;
solverOptions.CRmax = 0.3046 + 0.5189;
solverOptions.NPmin = 75;
solverOptions.crw = 0.2079;
solverOptions.fw = 0.3530;

% Common parameters
solverOptions.ftarget = 1e-8;
solverOptions.Display = 'off';
solverOptions.Noise = false;
solverOptions.EarlyStop = 'fitness';

% Parameter pairs to be analyzed
parameter_pairs = { ...
    'NPmin', ...
    'CRmin'; ...
    'NP', ...
    'CRmax'};

parameter_lbs = [...
    4,   ...    % NPmin
    eps; ...    % CRmin
    4, ...      % NP
    eps, ...    % CRmax
    ];

parameter_ubs = [...
    2000, ...   % NPmin
    1; ...      % CRmin
    2000, ...   % NP
    1, ...      % CRmax
    ];

% Sensitivity analysis options
sensitivity_options = cell(numel(parameter_pairs(1, :)), 1);
sensitivity_options{1}.XScale = 'log';
sensitivity_options{1}.YScale = 'log';
sensitivity_options{2}.XScale = 'default';
sensitivity_options{2}.YScale = 'default';

% Plot sensitivity analysis of paired parameters
for i = 1 : numel(parameter_pairs(1, :))    
    [SX, SY, SZ] = sensitivity_2D( ...
        solver, ...
        fitfun, ...
        lb, ...
        ub, ...
        maxfunevals, ...
        solverOptions, ...
        parameter_pairs(:, i), ...
        parameter_lbs(:, i), ...
        parameter_ubs(:, i), ...
        sensitivity_options{i});
    
    figure(i);    
    surf(SX, SY, SZ);
    set(gca, 'ZScale', 'log');
    
    xlabel(parameter_pairs{1, i});
    ylabel(parameter_pairs{2, i});
    zlabel('Solution Error');
    
    title(sprintf('Sensitivity Analysis of on %s', solver, fitfun), ...
        'Interpreter','none');
end
