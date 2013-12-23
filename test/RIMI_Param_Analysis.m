function RIMI_Param_Analysis(solver)
startTime = tic;

if nargin <= 0
	% 	solver = 'minmaxtcderand1binwoa';
	% 	solver = 'minmaxtcderand1bin';
	% 	solver = 'minmaxtcdebest1binwoa';
	% 	solver = 'minmaxtcdebest1bin';
	% 	solver = 'minmaxtcjadebinwoa';
	% 	solver = 'minmaxtcjadebin';
	% 	solver = 'minmaxdegl';
	solver = 'mmdeb1b_c';
end

if matlabpool('size') == 0
	matlabpool('open');
end

% Measure options
measureOptions.MaxFunEvals = 1e6;
measureOptions.Runs = 12;
measureOptions.FitnessFunctions = {'sainz_f1', 'sainz_f2', 'lu_f1'};
measureOptions.ConstraintFunctions = {'sainz_c1', 'sainz_c2', 'lu_c1'};
measureOptions.L1 = [-3.14, 0, -5];
measureOptions.L2 = [-3.14, 2, -5];
measureOptions.U1 = [3.14, 6, 5];
measureOptions.U2 = [3.14, 8, 5];

% Solver options
solverOptions1.dimensionFactor = 10;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.9;
solverOptions1.RecordPoint = 0;
solverOptions1.TolX = 0;
solverOptions1.TolFun = 0;
solverOptions1.TolStagnationIteration = 20;
solverOptions1.innerMaxIter = 200;
% solverOptions1.archiveSizeFactor = 4;
solverOptions2.dimensionFactor = 10;
solverOptions2.F = 0.7;
solverOptions2.CR = 0.9;
solverOptions2.TolStagnationIteration = 20;

reinitFactors = [0, 0.2, 0.4, 0.6, 0.8, 1];
migrateFactors = [0, 0.2, 0.4, 0.6, 0.8, 1];

for i = 1 : 6
	for j = 1 : 7 - i
		solverOptions1.reinitFactor = reinitFactors(i);
		solverOptions1.migrateFactor = migrateFactors(j);
		
		err = errors_cons_minmax(...
			solver, ...
			measureOptions, ...
			solverOptions1, ...
			solverOptions2);
		
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
		save(sprintf('err_con_minmax_%s_%s.mat', solver, datestr(now, 'yyyymmddHHMM')), ...
			'err', 'err_mean', 'err_std', 'sr', ...
			'measureOptions', 'solverOptions1', 'solverOptions2');
		fprintf('solver: %s\n', solver);
		fprintf('done\n');
	end
end
toc(startTime);
end
