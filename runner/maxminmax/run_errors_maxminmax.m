function run_errors_maxminmax(solver)
startTime = tic;

if nargin <= 0
	solver = 'mmmdeb1b_pce';
end

if matlabpool('size') <= 0
	matlabpool('open');
end

% Measure options
measureOptions.D1 = 1;
measureOptions.D2 = 1;
measureOptions.D3 = 1;
measureOptions.MaxFunEvals = 1e7;
measureOptions.Runs = 8;
measureOptions.L1 = -1 * ones(measureOptions.D1, 1);
measureOptions.L2 = -2 * ones(measureOptions.D2, 1);
measureOptions.L3 = -4 * ones(measureOptions.D3, 1);
measureOptions.U1 = 1 * ones(measureOptions.D1, 1);
measureOptions.U2 = 2 * ones(measureOptions.D2, 1);
measureOptions.U3 = 4 * ones(measureOptions.D3, 1);
measureOptions.FitnessFunctions = {...
	'maxminmax_f1', 'maxminmax_f8', 'maxminmax_f57', 'maxminmax_f64'};

% Solver options
solverOptions1.dimensionFactor = 10;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.9;
solverOptions1.RecordPoint = 0;
solverOptions1.TolX = 1e-9;
solverOptions1.TolFun = 0;
solverOptions1.TolStagnationIteration = 20;
solverOptions1.innerMaxIter = 200;
solverOptions1.migrateFactor = 0.7;
solverOptions1.zeta = 1e-4;
solverOptions2.dimensionFactor = 10;
solverOptions2.F = 0.7;
solverOptions2.CR = 0.9;
solverOptions3.dimensionFactor = 10;
solverOptions3.F = 0.7;
solverOptions3.CR = 0.9;

err = errors_maxminmax(solver, measureOptions, solverOptions1, ...
	solverOptions2, solverOptions3);

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
save(sprintf('err_%s_%s.mat', solver, datestr(now, 'yyyymmddHHMM')), ...
	'err', 'err_mean', 'err_std', 'sr', ...
	'measureOptions', 'solverOptions1', 'solverOptions2', 'solverOptions3');
fprintf('solver: %s\n', solver);
fprintf('done\n');
toc(startTime);
end
