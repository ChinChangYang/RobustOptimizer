function ret = errors_cec2013(solver, measureOptions, solverOptions)
% ERRORS Errors of a solver over a set of test functions

% Handle input options
if nargin <= 2
	solverOptions = [];
end

% Handle measure options
if nargin <= 1
	measureOptions = [];
end

if ~isfield(measureOptions, 'Dimension')
	D = 10;
else
	D = measureOptions.Dimension;
end

if ~isfield(measureOptions, 'Runs')
	nRuns = 51;
else
	nRuns = measureOptions.Runs;
end

if ~isfield(measureOptions, 'FitnessFunctions')
	fitfuns = {'cec13_f1', 'cec13_f2', 'cec13_f3', 'cec13_f4', ...
		'cec13_f5', 'cec13_f6', 'cec13_f7', 'cec13_f8', 'cec13_f9', ...
		'cec13_f10', 'cec13_f11', 'cec13_f12', 'cec13_f13', ...
		'cec13_f14', 'cec13_f15', 'cec13_f16', 'cec13_f17', ...
		'cec13_f18', 'cec13_f19', 'cec13_f20', 'cec13_f21', ...
		'cec13_f22', 'cec13_f23', 'cec13_f24', 'cec13_f25', ...
		'cec13_f26', 'cec13_f27', 'cec13_f28'};
else
	fitfuns = measureOptions.FitnessFunctions;
end

if ~isfield(measureOptions, 'MaxFunEvalSet')
	maxfunevalset = round(D .* 1e4);
else
	maxfunevalset = measureOptions.MaxFunEvalSet;
end

if ~isfield(measureOptions, 'LowerBounds')
	lb = -100 * ones(D, 1);
else
	lb = eval(measureOptions.LowerBounds);
end

if ~isfield(measureOptions, 'UpperBounds')
	ub = 100 * ones(D, 1);
else
	ub = eval(measureOptions.UpperBounds);
end

% Statistic variables
allfmin = zeros(length(maxfunevalset), length(fitfuns), nRuns);

% Main optimizer
rand('state', sum(100*clock)); %#ok<RAND>
for iFitfun = 1 : length(fitfuns)
	fitfun = fitfuns{iFitfun};
	for iMaxfuneval = 1 : length(maxfunevalset)
		% Reset default random stream
		maxfunevals = maxfunevalset(iMaxfuneval);
        parfor iRuns = 1 : nRuns
			fprintf('fitfun: %s, maxfunevals: %d, runs: %d\n', ...
				fitfun, maxfunevals, iRuns);
            [~, allfmin(iMaxfuneval, iFitfun, iRuns)] = ...
                feval(solver, fitfun, lb, ub, maxfunevals, ...
				solverOptions);
        end
	end
end

% Mean solution errors
% meanErrs = mean(allfmin, 3);

% Return results
ret = allfmin;
end
