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

% Sensitivity analysis parameters
parameter_set = { ...
    'NP', ...
    'F', ...
    'CR', ...
    'ER', ...
    'p', ...
    'H', ...
    'Q', ...
    'Ar', ...
    'cw', ...
    'erw', ...
    'CRmin', ...
    'CRmax', ...
    'NPmin', ...
    'crw', ...
    'fw', ...
    };

parameter_lbs = [...
    solverOptions.NPmin, ... % NP
    eps, ... % F
    eps, ... % CR
    eps, ... % ER
    eps, ... % p
    2,   ... % H
    0,   ... % Q
    0,   ... % Ar
    eps, ... % cw
    eps, ... % erw
    eps, ... % CRmin
    solverOptions.CRmin, ... % CRmax
    4,   ... % NPmin
    eps, ... % crw
    eps, ... % fw
    ];

parameter_ubs = [...
    3000, ... % NP
    1, ... % F
    1, ... % CR
    1, ... % ER
    1, ... % p
    3000, ... % H
    3000, ... % Q
    10,   ... % Ar
    1, ... % cw
    1, ... % erw
    solverOptions.CRmax, ... % CRmin
    1, ... % CRmax
    solverOptions.NP, ... % NPmin
    1, ... % crw
    1, ... % fw
    ];

% Estimate sensitivity
[SX, SY, meanfmin, stdfmin] = sensitivity( ...
    solver, ...
    fitfun, ...
    lb, ...
    ub, ...
    maxfunevals, ...
    solverOptions, ...
    parameter_set, ...
    parameter_lbs, ...
    parameter_ubs);

% Plot sensitivity analysis results
for i = 1 : numel(parameter_set)
    figure(i);
        
    plot(SX(:, i), SY(:, i));
    hold on;
    plot(SX(:, i), repmat(meanfmin - stdfmin, 1, 5), 'r');
    plot(SX(:, i), repmat(meanfmin + stdfmin, 1, 5), 'r');
    hold off;
    
    title(sprintf('Sensitivity Analysis of %s on %s', solver, fitfun), ...
        'Interpreter','none');
    
    xlabel(parameter_set{i});
    ylabel('Solution Error');
end