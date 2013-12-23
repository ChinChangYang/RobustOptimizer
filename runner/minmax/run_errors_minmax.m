function run_errors_minmax(solver)
startTime = tic;

if nargin <= 0
	solver = 'mmjade_pce';
end

if matlabpool('size') == 0
    matlabpool('open');
end

% Measure options
measureOptions.D1 = 2;
measureOptions.D2 = 2;
measureOptions.MaxFunEvals = 2 * 2 * 1e6;
measureOptions.Runs = 12;
measureOptions.L1 = -2 * ones(measureOptions.D1, 1);
measureOptions.L2 = -4 * ones(measureOptions.D2, 1);
measureOptions.U1 = 2 * ones(measureOptions.D1, 1);
measureOptions.U2 = 4 * ones(measureOptions.D2, 1);
measureOptions.FitnessFunctions = {...
	'fminmax_f1', 'fminmax_f2', 'fminmax_f3', ...
	'fminmax_f4', 'fminmax_f5', 'fminmax_f6', ...
	'fminmax_f7', 'fminmax_f8'};

% Solver options
solverOptions1.dimensionFactor = 15;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.9;
solverOptions1.RecordPoint = 0;
solverOptions1.TolX = 0;
solverOptions1.TolFun = 0;
solverOptions1.TolStagnationIteration = 20;
solverOptions1.innerMaxIter = 200;
solverOptions1.migrateFactor = 0.7;
solverOptions2.dimensionFactor = 15;
solverOptions2.F = 0.7;
solverOptions2.CR = 0.9;
solverOptions2.TolStagnationIteration = 20;

err = errors_minmax(solver, measureOptions, solverOptions1, solverOptions2);
[~, ~, nRuns] = size(err);
err_mean = mean(err, 3);
err_std = std(err, [], 3);
sr = sum(err <= 0, 3) ./ nRuns;
fprintf('err = \n');
disp(err);
fprintf('err_mean = \n');
disp(err_mean);
fprintf('err_std = \n');
disp(err_std);
fprintf('sr = \n');
disp(sr);
save(sprintf('err_minmax_%s_%s.mat', solver, datestr(now, 'yyyymmddHHMM')), ...
	'err', 'err_mean', 'err_std', 'sr', ...
	'measureOptions', 'solverOptions1', 'solverOptions2');
fprintf('solver: %s\n', solver);
fprintf('done\n');
toc(startTime);
end
