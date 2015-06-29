% Common control
clear;
close all;

% Experiment settings
solver = 'SPS_L_SHADE_EIG';
fitfun = 'cec15_f3';
D = 10;
maxfunevals = D * 1e4;
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);

% Parameter baseline
solverOptions.NP = 9 + 155;
solverOptions.F = 0.4085;
solverOptions.CR = 0.1980;
solverOptions.ER = 0.3172;
solverOptions.p = 0.2230;
solverOptions.H = 431;
solverOptions.Q = 151;
solverOptions.Ar = 1.7238;
solverOptions.cw = 0.1720;
solverOptions.erw = 0.1717;
solverOptions.CRmin = 0.0100;
solverOptions.CRmax = 0.0100 + 0.9248;
solverOptions.NPmin = 9;
solverOptions.crw = 0.9933;
solverOptions.fw = 0.4899;

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
    1,   ... % H
    0,   ... % Q
    1,   ... % Ar
    eps, ... % cw
    eps, ... % erw
    eps, ... % CRmin
    solverOptions.CRmin, ... % CRmax
    4,   ... % NPmin
    eps, ... % crw
    eps, ... % fw
    ];

parameter_ubs = [...
    2500, ... % NP
    1, ... % F
    1, ... % CR
    1, ... % ER
    1, ... % p
    2500, ... % H
    2500, ... % Q
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
[SX, SY, SMF, SSF] = sensitivity( ...
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
for i = 1 : numel(SMF)
    figure(i);
        
    plot(SX(:, i), SY(:, i));
    hold on;
    plot(SX(:, i), repmat(SMF(i) - SSF(i), 1, 5), 'r');
    plot(SX(:, i), repmat(SMF(i) + SSF(i), 1, 5), 'r');
    hold off;
    
    title(sprintf('Sensitivity Analysis of %s on %s', solver, fitfun), ...
        'Interpreter','none');
    
    xlabel(parameter_set{i});
    ylabel('Solution Error');
end