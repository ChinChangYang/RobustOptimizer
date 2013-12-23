function ret = computeerrors(errorMeasure, solver, genotype, measureOptions)
%COMPUTEERRORS Compute errors
if nargin <= 3
    measureOptions = [];
end

defaultOptions.D = 10;
defaultOptions.Runs = 2;
defaultOptions.MaxFunEvals = 1e5;
defaultOptions.LowerBounds = -100;
defaultOptions.UpperBounds = 100;
defaultOptions.RecordPoint = 0;

optionNames = fieldnames(defaultOptions);

for i = 1 : numel(optionNames)
    if ~isfield(measureOptions, optionNames{i})
        defaultOption = getfield(defaultOptions, optionNames{i}); %#ok<GFLD>
        measureOptions = setfield(measureOptions, optionNames{i}, defaultOption); %#ok<SFLD>
    end
end

D = measureOptions.D;
nRuns = measureOptions.Runs;
maxfunevals = measureOptions.MaxFunEvals;
if isscalar(measureOptions.LowerBounds)
    lb = measureOptions.LowerBounds * ones(D, 1);
else
    lb = measureOptions.LowerBounds;
end
if isscalar(measureOptions.UpperBounds)
    ub = measureOptions.UpperBounds * ones(D, 1);
else
    ub = measureOptions.UpperBounds;
end

solverOptions = genotype2options(genotype, solver);

% Statistic variables
allfmin = zeros(1, nRuns);

% Main optimizer
for iRuns = 1 : nRuns
	[~, allfmin(iRuns)] = ...
		feval(solver, errorMeasure, lb, ub, maxfunevals, ...
		solverOptions);
end

% Mean of the minimal functions
ret = mean(allfmin);
end

