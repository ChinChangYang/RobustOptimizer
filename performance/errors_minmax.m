function ret = errors_minmax(solver, measureOptions, solverOptions1, ...
	solverOptions2)
%ERRORS_MINMAX Errors of a solver within a set of min-max test functions
if nargin <= 0
	solver = 'minmaxtcjadebin';
end

% Handle measure options
if nargin <= 1
	measureOptions = [];
end

% Handle input options
if nargin <= 2
	solverOptions1 = [];
end

if nargin <= 3
	solverOptions2 = [];
end

% Default measure options
defaultMeasureOptions.FitnessFunctions = {...
	'fminmax_f1', 'fminmax_f2', 'fminmax_f3', ...
	'fminmax_f4', 'fminmax_f5', 'fminmax_f6', ...
	'fminmax_f7', 'fminmax_f8'};
defaultMeasureOptions.D1 = 2; % Dimension of minimizers
defaultMeasureOptions.D2 = 2; % Dimension of maximizers
defaultMeasureOptions.Runs = 30;
defaultMeasureOptions.MaxFunEvals = 2e4;
defaultMeasureOptions.L1 = -2 * ones(defaultMeasureOptions.D1, 1);
defaultMeasureOptions.L2 = -4 * ones(defaultMeasureOptions.D2, 1);
defaultMeasureOptions.U1 = 2 * ones(defaultMeasureOptions.D1, 1);
defaultMeasureOptions.U2 = 4 * ones(defaultMeasureOptions.D2, 1);
measureOptions = setdefoptions(measureOptions, defaultMeasureOptions);

% Default solver options
defaultSolverOptions.Display = 'off';
defaultSolverOptions.RecordPoint = 0;
solverOptions1 = setdefoptions(solverOptions1, defaultSolverOptions);
solverOptions2 = setdefoptions(solverOptions2, defaultSolverOptions);

% Bounds
lb1 = measureOptions.L1;
lb2 = measureOptions.L2;
ub1 = measureOptions.U1;
ub2 = measureOptions.U2;

% Initialize experimental results
results = zeros(numel(measureOptions.FitnessFunctions), ...
	numel(measureOptions.MaxFunEvals), ...
	measureOptions.Runs);

for iFitfun = 1 : numel(measureOptions.FitnessFunctions)
	fitfun = measureOptions.FitnessFunctions{iFitfun};
	for iMaxfunevals = 1 : numel(measureOptions.MaxFunEvals)
		maxfunevals = measureOptions.MaxFunEvals(iMaxfunevals);
		parfor iRuns = 1 : measureOptions.Runs
			rng(iRuns, 'twister');
			
			[xminmax, ~, ~] = ...
				feval(solver, fitfun, maxfunevals, ...
				lb1, ub1, lb2, ub2, solverOptions1, solverOptions2);
			
			results(iFitfun, iMaxfunevals, iRuns) = ...
				norm(xminmax - feval(fitfun));
			
			if results(iFitfun, iMaxfunevals, iRuns) < 1e-6
				results(iFitfun, iMaxfunevals, iRuns) = 0;
			end
			
			fprintf('fitfun: %s, maxfunevals: %d, runs: %d, error: %.4E\n', ...
				fitfun, maxfunevals, iRuns, results(iFitfun, iMaxfunevals, iRuns));
		end
	end
end

ret = results;
end
