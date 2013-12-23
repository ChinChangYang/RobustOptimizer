function ret = errors_bbob12(solver, measureOptions, solverOptions)
% ERRORS Errors of a solver over a set of test functions

% Handle input options
if nargin <= 2
	solverOptions = [];
end

% Handle measure options
if nargin <= 1
	measureOptions = [];
end

defMeasureOptions.Dimension = 10;
defMeasureOptions.Runs = 51;
defMeasureOptions.MaxFunEvalSet = 'round(D * 1e4)';
defMeasureOptions.LowerBounds = '-5 * ones(D, 1)';
defMeasureOptions.UpperBounds = '5 * ones(D, 1)';
defMeasureOptions.FitnessFunctions = {...
	'bbob12_f1', 'bbob12_f2', 'bbob12_f3', 'bbob12_f4', ...
	'bbob12_f5', 'bbob12_f6', 'bbob12_f7', 'bbob12_f8', 'bbob12_f9', ...
	'bbob12_f10', 'bbob12_f11', 'bbob12_f12', 'bbob12_f13', ...
	'bbob12_f14', 'bbob12_f15', 'bbob12_f16', 'bbob12_f17', ...
	'bbob12_f18', 'bbob12_f19', 'bbob12_f20', 'bbob12_f21', ...
	'bbob12_f22', 'bbob12_f23', 'bbob12_f24'};
measureOptions = setdefoptions(measureOptions, defMeasureOptions);

D = measureOptions.Dimension; %#ok<NASGU>
nRuns = measureOptions.Runs;
fitfuns = measureOptions.FitnessFunctions;
if ischar(measureOptions.MaxFunEvalSet)
	maxfunevalset = eval(measureOptions.MaxFunEvalSet);
else
	maxfunevalset = measureOptions.MaxFunEvalSet;
end
lb = eval(measureOptions.LowerBounds);
ub = eval(measureOptions.UpperBounds);

% Statistic variables
allfmin = zeros(length(maxfunevalset), length(fitfuns), nRuns);

% Main optimizer
rand('state', sum(100*clock)); %#ok<RAND>
for iFitfun = 1 : length(fitfuns)
	fitfun = fitfuns{iFitfun};
	for iMaxfuneval = 1 : length(maxfunevalset)
		maxfunevals = maxfunevalset(iMaxfuneval);
		parfor iRuns = 1 : nRuns
			fprintf('fitfun: %s, maxfunevals: %d, runs: %d\n', ...
				fitfun, maxfunevals, iRuns);
			[~, allfmin(iMaxfuneval, iFitfun, iRuns)] = ...
				feval(solver, fitfun, lb, ub, maxfunevals, ...
				solverOptions);
		end
		
		fopt = feval(fitfun, 'xopt');
		allfmin(iMaxfuneval, iFitfun, :) = ...
			allfmin(iMaxfuneval, iFitfun, :) - fopt;		
	end
end

% Tolerant errors
for iFitfun = 1 : length(fitfuns)
	for iMaxfuneval = 1 : length(maxfunevalset)
		for iRuns = 1 : nRuns
			if allfmin(iMaxfuneval, iFitfun, iRuns) < 1e-8
				allfmin(iMaxfuneval, iFitfun, iRuns) = 0;
			end
		end
	end
end

% Return results

ret = allfmin;
end
