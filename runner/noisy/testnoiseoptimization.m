function testnoiseoptimization
%TESTNOISEOPTIMIZATION Test noise optimization
addprojectpath;
solver = 'projadeeig';
noisyfitfun = 'ellipsoidalrotnoise';
noisefreefitfun = 'ellipsoidalrot';
optimfval = 0;
D = 5;
maxfunevals = 1e5;
solverOptions.dimensionFactor = 6.3302e+00;
solverOptions.TolX = 0;
solverOptions.TolFun = 0;
solverOptions.Restart = 0;
solverOptions.Display = 'off';
solverOptions.Noise = true;
solverOptions.StopOnStagnation = 'off';
solverOptions.StopOnEqualFunctionValues = false;
solverOptions.StopOnWarnings = 'no';
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);
fmin = zeros(1, 10);
for i = 1 : numel(fmin)
	[xmin, ~, ~] = ...
		feval(solver, noisyfitfun, lb, ub, maxfunevals, solverOptions);
	fmin(i) = feval(noisefreefitfun, xmin);
end

fprintf('Solver: %s\n', solver);
fprintf('Maximal function evaluations: %.1E\n', maxfunevals);
fprintf('Mean of solution errors: %.4E\n', mean(abs(fmin - optimfval)));
fprintf('Std. of solution errors: %.4E\n', std(abs(fmin - optimfval)));
end
