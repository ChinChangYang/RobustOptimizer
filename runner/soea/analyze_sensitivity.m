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
% parameter = 'NP';
% parameter_lb = solverOptions.NPmin;
% parameter_ub = 2500;

% parameter = 'NPmin';
% parameter_lb = 4;
% parameter_ub = solverOptions.NP;

% parameter = 'Ar';
% parameter_lb = 1;
% parameter_ub = 4;

parameter = 'erw';
parameter_lb = eps;
parameter_ub = 1;

if parameter_ub < solverOptions.(parameter)
    parameter_ub = max(0, solverOptions.(parameter) * 2);
end

% Estimate sensitivity
[X, Y, meanfmin, stdfmin] = sensitivity( ...
    solver, ...
    fitfun, ...
    lb, ...
    ub, ...
    maxfunevals, ...
    solverOptions, ...
    parameter, ...
    parameter_lb, ...
    parameter_ub);

% Plot sensitivity analysis result
plot(X, Y);
hold on;
plot(X, repmat(meanfmin - stdfmin, 1, 5), 'r');
plot(X, repmat(meanfmin + stdfmin, 1, 5), 'r');
hold off;

title(sprintf('Sensitivity Analysis of %s on %s', solver, fitfun), ...
    'Interpreter','none');

xlabel(parameter);
ylabel('Solution Error');