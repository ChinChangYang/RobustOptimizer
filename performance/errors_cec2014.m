function ret = errors_cec2014(solver, measureOptions, solverOptions)
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
	fitfuns = {'cec14_f1', 'cec14_f2', 'cec14_f3', 'cec14_f4', ...
		'cec14_f5', 'cec14_f6', 'cec14_f7', 'cec14_f8', 'cec14_f9', ...
		'cec14_f10', 'cec14_f11', 'cec14_f12', 'cec14_f13', ...
		'cec14_f14', 'cec14_f15', 'cec14_f16', 'cec14_f17', ...
		'cec14_f18', 'cec14_f19', 'cec14_f20', 'cec14_f21', ...
		'cec14_f22', 'cec14_f23', 'cec14_f24', 'cec14_f25', ...
		'cec14_f26', 'cec14_f27', 'cec14_f28', 'cec14_f29', ...
		'cec14_f30'};
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
