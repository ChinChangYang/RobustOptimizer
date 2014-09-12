function ret = errors_cons_minmax(solver, measureOptions, solverOptions1, ...
	solverOptions2)
%ERRORS_CONS_MINMAX Errors of a solver within a set of min-max test
%functions with constraints
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
	'sainz_f1', 'sainz_f2', 'lu_f1'};
defaultMeasureOptions.ConstraintFunctions = {...
	'sainz_c1', 'sainz_c2', 'lu_c1'};
defaultMeasureOptions.D1 = 1; % Dimension of minimizers
defaultMeasureOptions.D2 = 1; % Dimension of maximizers
defaultMeasureOptions.Runs = 30;
defaultMeasureOptions.MaxFunEvals = 2e4;
defaultMeasureOptions.L1 = [-3.14, 0, -5];
defaultMeasureOptions.L2 = [-3.14, 2, -5];
defaultMeasureOptions.U1 = [3.14, 6, 5];
defaultMeasureOptions.U2 = [3.14, 8, 5];
measureOptions = setdefoptions(measureOptions, defaultMeasureOptions);

% Default solver options
defaultSolverOptions.Display = 'off';
defaultSolverOptions.RecordPoint = 0;
solverOptions1 = setdefoptions(solverOptions1, defaultSolverOptions);
solverOptions2 = setdefoptions(solverOptions2, defaultSolverOptions);

% Initialize experimental results
results = zeros(numel(measureOptions.FitnessFunctions), ...
	numel(measureOptions.MaxFunEvals), ...
	measureOptions.Runs);

for iFitfun = 1 : numel(measureOptions.FitnessFunctions)
	fitfun = measureOptions.FitnessFunctions{iFitfun};
	nonlcon = measureOptions.ConstraintFunctions{iFitfun};
	solverOptions1.nonlcon = nonlcon;
	
	% Bounds
	lb1 = measureOptions.L1(:, iFitfun) * ones(measureOptions.D1, 1);
	lb2 = measureOptions.L2(:, iFitfun) * ones(measureOptions.D2, 1);
	ub1 = measureOptions.U1(:, iFitfun) * ones(measureOptions.D1, 1);
	ub2 = measureOptions.U2(:, iFitfun) * ones(measureOptions.D2, 1);

	for iMaxfunevals = 1 : numel(measureOptions.MaxFunEvals)
		maxfunevals = measureOptions.MaxFunEvals(iMaxfunevals);
		parfor iRuns = 1 : measureOptions.Runs
			rng(iRuns, 'twister');
			
			[xbest1, xbest2, ~] = ...
				feval(solver, fitfun, maxfunevals, ...
				lb1, ub1, lb2, ub2, solverOptions1, solverOptions2);
			
			xbest = [xbest1; xbest2];
			optima = feval(fitfun);
			[~, n_optima] = size(optima);
			distance = zeros(n_optima, 1);
			
			for i_optima = 1 : n_optima
				distance(i_optima) = norm(xbest - optima(:, i_optima));
			end
			
			[min_distance, min_index] = min(distance);			
			results(iFitfun, iMaxfunevals, iRuns) = min_distance;
			
			if all(abs(xbest - optima(:, min_index)) < 1e-6)
				results(iFitfun, iMaxfunevals, iRuns) = 0;
			end
			
			fprintf('fitfun: %s, maxfunevals: %d, runs: %d, error: %.4E\n', ...
				fitfun, maxfunevals, iRuns, results(iFitfun, iMaxfunevals, iRuns));
		end
	end
end

ret = results;
end
