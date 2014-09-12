function run_errors_cons_minmax2(solver)
startTime = tic;

if nargin <= 0
	solver = 'mmjade_pce';
end

if matlabpool('size') == 0
    matlabpool('open');
end

% Measure options
measureOptions.MaxFunEvals = 2e7;
measureOptions.Runs = 8;
measureOptions.FitnessFunctions = {...
	'fminmax_f1', ...
	'fminmax_f2'};
measureOptions.ConstraintFunctions = {...
	'fminmax_c1', ...
	'fminmax_c2'};
measureOptions.D1 = 2;
measureOptions.D2 = 2;
measureOptions.L1 = [-2, -2];
measureOptions.L2 = [-4, -4];
measureOptions.U1 = [2, 2];
measureOptions.U2 = [4, 4];

% Solver options
solverOptions1.NP = 30;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.9;
solverOptions1.RecordPoint = 0;
solverOptions1.TolStagnationIteration = 20;
solverOptions1.innerMaxIter = 200;
solverOptions1.migrateFactor = 0.7;
solverOptions1.ConstraintHandling = 'EpsilonMethod';
solverOptions1.EpsilonValue = 1e-6;
solverOptions1.EarlyStop = 'auto';
solverOptions2.NP = 30;
solverOptions2.F = 0.7;
solverOptions2.CR = 0.9;
solverOptions2.TolStagnationIteration = 20;
solverOptions2.ConstraintHandling = 'EpsilonMethod';
solverOptions2.EpsilonValue = 1e-6;
solverOptions2.EarlyStop = 'auto';

err = errors_cons_minmax(solver, measureOptions, solverOptions1, solverOptions2);
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
save(sprintf('err_con_minmax2_%s_%s.mat', solver, datestr(now, 'yyyymmddHHMM')), ...
	'err', 'err_mean', 'err_std', 'sr', ...
	'measureOptions', 'solverOptions1', 'solverOptions2');
fprintf('solver: %s\n', solver);
fprintf('done\n');
toc(startTime);
end
