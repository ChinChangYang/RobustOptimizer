%RUN_BUILTIN_PSO Run built-in PSO global optimization
clear;
close all;
startTime = tic;

% rng('default');
solver = 'particleswarm';
fitfun = @(x) cec15_f1(x);
D = 10;
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);
% lb = -6.4 * ones(D, 1);
% ub = 6.35 * ones(D, 1);
% lb = zeros(20, 1);
% ub = 4 * pi * ones(20, 1);
% [lb, ub] = getlimit_messenger;
% [lb, ub] = getlimit_cassini2;
hybridopts = optimoptions('fmincon', 'FunctionTolerance', 1e-8);
solverOptions = optimoptions('particleswarm', ...
    'FunctionTolerance', 1e-8, ...
    'HybridFcn', {@fmincon, hybridopts});

[xmin, fmin, exitflag, output] = ...
    feval(solver, fitfun, D, lb, ub, solverOptions);

fprintf('xmin =\n');
disp(xmin);
fprintf('fmin =\n');
fprintf('    %E\n', fmin);
fprintf('funccount = %d\n', output.funccount);

toc(startTime);
