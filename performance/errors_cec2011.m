function ret = errors_cec2011(solver, measureOptions, solverOptions)
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
	D = [6, 20, 12, 26, 22];
else
	D = measureOptions.Dimension;
end

if ~isfield(measureOptions, 'Runs')
	nRuns = 51;
else
	nRuns = measureOptions.Runs;
end

if ~isfield(measureOptions, 'FitnessFunctions')
	fitfuns = {'cec11_f1', 'cec11_f7', 'cec11_f10', 'cec11_f12', ...
		'cec11_f13'};
else
	fitfuns = measureOptions.FitnessFunctions;
end

if ~isfield(measureOptions, 'MaxFunEvalSet')
	maxfunevalset = round(D .* 1e4);
else
	maxfunevalset = measureOptions.MaxFunEvalSet;
end

% Statistic variables
allfmin = zeros(length(fitfuns), nRuns);

% Main optimizer
rand('state', sum(100*clock)); %#ok<RAND>
for iFitfun = 1 : length(fitfuns)
	maxfunevals = maxfunevalset(iFitfun);
	fitfun = fitfuns(iFitfun);
	fitfun = fitfun{:};
	if strcmp(fitfun, 'cec11_f1')
		lb = -6.4 * ones(6, 1);
		ub = 6.35 * ones(6, 1);
	elseif strcmp(fitfun, 'cec11_f7')
		lb = zeros(20, 1);
		ub = 4 * pi * ones(20, 1);
	elseif strcmp(fitfun, 'cec11_f10')
		lb = [0.2 * ones(6, 1); -180 * ones(6, 1)];
		ub = [ones(6, 1); 180 * ones(6, 1)];
	elseif strcmp(fitfun, 'cec11_f12')
		[lb, ub] = getlimit_messenger;
	elseif strcmp(fitfun, 'cec11_f13')
		[lb, ub] = getlimit_cassini2;
	else
		fprintf('Exception: Unknown fitfun: %s\n', fitfun);
		fprintf('Aborted!\n');
		continue;
	end
	
	% Reset default random stream
	for iRuns = 1 : nRuns
		
		fprintf('fitfun: %s, maxfunevals: %d, runs: %d\n', ...
			fitfun, maxfunevals, iRuns);
		[~, allfmin(iFitfun, iRuns)] = ...
			feval(solver, fitfun, lb, ub, maxfunevals, ...
			solverOptions);
	end
end

% Mean solution errors
% meanErrs = mean(allfmin, 3);

% Return results
ret = allfmin;
end
