function ret = errors_maxminmax(solver, measureOptions, solverOptions1, ...
	solverOptions2, solverOptions3)
%ERRORS_MINMAX Errors of a solver within a set of min-max test functions

if nargin <= 0
	solver = 'maxminmaxtcjade';
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

if nargin <= 4
	solverOptions3 = [];
end

% Reset random number generator
rng(0, 'twister');

% Default measure options
defaultMeasureOptions.FitnessFunctions = {...
	'maxminmax_f1', 'maxminmax_f8', 'maxminmax_f57', 'maxminmax_f64'};

defaultMeasureOptions.D1 = 2;
defaultMeasureOptions.D2 = 2;
defaultMeasureOptions.D3 = 2;
defaultMeasureOptions.Runs = 30;
defaultMeasureOptions.MaxFunEvals = 2e4;
defaultMeasureOptions.L1 = -1 * ones(defaultMeasureOptions.D1, 1);
defaultMeasureOptions.L2 = -2 * ones(defaultMeasureOptions.D2, 1);
defaultMeasureOptions.L3 = -4 * ones(defaultMeasureOptions.D3, 1);
defaultMeasureOptions.U1 = 1 * ones(defaultMeasureOptions.D1, 1);
defaultMeasureOptions.U2 = 2 * ones(defaultMeasureOptions.D2, 1);
defaultMeasureOptions.U3 = 4 * ones(defaultMeasureOptions.D3, 1);
measureOptions = setdefoptions(measureOptions, defaultMeasureOptions);

% Default solver options
defaultSolverOptions.Display = 'off';
defaultSolverOptions.RecordPoint = 0;
solverOptions1 = setdefoptions(solverOptions1, defaultSolverOptions);
solverOptions2 = setdefoptions(solverOptions2, defaultSolverOptions);
solverOptions3 = setdefoptions(solverOptions3, defaultSolverOptions);

% Bounds
lb1 = measureOptions.L1;
ub1 = measureOptions.U1;
lb2 = measureOptions.L2;
ub2 = measureOptions.U2;
lb3 = measureOptions.L3;
ub3 = measureOptions.U3;

% Initialize experimental results
results = zeros(numel(measureOptions.FitnessFunctions), ...
	numel(measureOptions.MaxFunEvals), ...
	measureOptions.Runs);

for iFitfun = 1 : numel(measureOptions.FitnessFunctions)
	fitfun = measureOptions.FitnessFunctions{iFitfun};
	shift = feval(fitfun);
	
	for iMaxfunevals = 1 : numel(measureOptions.MaxFunEvals)
		
		maxfunevals = measureOptions.MaxFunEvals(iMaxfunevals);
		
		parfor iRuns = 1 : measureOptions.Runs	
			rng(iRuns, 'twister');		
			
			[xminmax, ~, ~] = ...
				feval(solver, fitfun, maxfunevals, ...
				lb1, ub1, lb2, ub2, lb3, ub3, ...
				solverOptions1, solverOptions2, solverOptions3);
			
			results(iFitfun, iMaxfunevals, iRuns) = norm(xminmax - shift);
			
			if results(iFitfun, iMaxfunevals, iRuns) <= 1e-6
				results(iFitfun, iMaxfunevals, iRuns) = 0;
			end
			
			fprintf('fitfun: %s, maxfunevals: %d, runs: %d, error: %.4E\n', ...
				fitfun, maxfunevals, iRuns, results(iFitfun, iMaxfunevals, iRuns));
		end
	end
end

ret = results;
end
