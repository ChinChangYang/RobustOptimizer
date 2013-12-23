function run_noisy_experiment
%RUN_NOISY_EXPERIMENT Run noisy experiment
addprojectpath;
solvers = {'projadeeig', 'debest1eig', 'jadeeig', 'sadeeig'};
maxfunevals = {1e5, 1e6};
noisyfitfun = 'ellipsoidalrotnoise';
noisefreefitfun = 'ellipsoidalrot';
optimfval = 0;
D = 5;
solverOptions.TolX = 0;
solverOptions.TolFun = 0;
solverOptions.Restart = 0;
solverOptions.Display = 'off';
solverOptions.Noise = true;
solverOptions.SampleFactor = 1.01;
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);
for iSolver = 1 : numel(solvers)
	for iMaxFunEvals = 1 : numel(maxfunevals)
		solver = solvers{iSolver};
		maxfuneval = maxfunevals{iMaxFunEvals};
		fmin = zeros(1, 10);
		for i = 1 : numel(fmin)
			[xmin, ~, ~] = ...
				feval(solver, noisyfitfun, lb, ub, maxfuneval, solverOptions);
			fmin(i) = feval(noisefreefitfun, xmin);
		end
		
		fprintf('Solver: %s\n', solver);
		fprintf('Maximal function evaluations: %.4E\n', maxfuneval);
		fprintf('Mean of solution errors: %.4E\n', mean(abs(fmin - optimfval)));
		fprintf('Std. of solution errors: %.4E\n', std(abs(fmin - optimfval)));
	end
end
end

