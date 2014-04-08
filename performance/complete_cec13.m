function [allout, allfvals, allfes, T0, T1, T2] = complete_cec13(...
	solver, ...
	measureOptions, ...
	solverOptions)
% COMPLETE_CEC13 complete CEC'13 experiments

% Deal with input arguments
if nargin <= 0
	solver = 'derand1bin';
end

if nargin <= 1
	measureOptions = [];
end

if nargin <= 2
	solverOptions = [];
end

defaultMeasureOptions.Dimension = 30;
defaultMeasureOptions.Runs = 51;
defaultMeasureOptions.FitnessFunctions = ...
	{'cec13_f1', 'cec13_f2', 'cec13_f3', 'cec13_f4', ...
	'cec13_f5', 'cec13_f6', 'cec13_f7', 'cec13_f8', 'cec13_f9', ...
	'cec13_f10', 'cec13_f11', 'cec13_f12', 'cec13_f13', ...
	'cec13_f14', 'cec13_f15', 'cec13_f16', 'cec13_f17', ...
	'cec13_f18', 'cec13_f19', 'cec13_f20', 'cec13_f21', ...
	'cec13_f22', 'cec13_f23', 'cec13_f24', 'cec13_f25', ...
	'cec13_f26', 'cec13_f27', 'cec13_f28'};
defaultMeasureOptions.MaxFunEvals = 3e5;
defaultMeasureOptions.LowerBounds = -100;
defaultMeasureOptions.UpperBounds = 100;
measureOptions = setdefoptions(measureOptions, defaultMeasureOptions);
D		= measureOptions.Dimension;
runs	= measureOptions.Runs;
fitfuns = measureOptions.FitnessFunctions;
maxfes	= measureOptions.MaxFunEvals;
lb		= measureOptions.LowerBounds * ones(D, 1);
ub		= measureOptions.UpperBounds * ones(D, 1);

defaultSolverOptions.dimensionFactor = 5;
defaultSolverOptions.NP = defaultSolverOptions.dimensionFactor * D;
defaultSolverOptions.RecordPoint = 21;
defaultSolverOptions.ftarget = 1e-8;
defaultSolverOptions.TolX = 0;
defaultSolverOptions.TolStagnationIteration = Inf;
defaultSolverOptions.initial.X = ...
	repmat(lb, 1, defaultSolverOptions.NP) + ...
	repmat(ub - lb, 1, defaultSolverOptions.NP) .* ...
	lhsdesign(defaultSolverOptions.NP, D, 'iteration', 100)';

solverOptions = setdefoptions(solverOptions, defaultSolverOptions);
progress = solverOptions.RecordPoint;

% Statistic variables
allfes = zeros(progress, runs, numel(fitfuns));
allfvals = zeros(progress, runs, numel(fitfuns));
allout = cell(runs, numel(fitfuns));

% Collect function values and output
pause(rand); % for the reseeding in the next line
rand('state', sum(100*clock)); %#ok<RAND>
for ifitfun = 1 : length(fitfuns)
	fitfun = fitfuns(ifitfun);
	fitfun = fitfun{:};
	parfor iruns = 1 : runs
		fprintf('fitfun: %s, maxfes: %d, iruns: %d\n', ...
			fitfun, maxfes, iruns);
		[~, ~, out] = ...
			feval(solver, fitfun, lb, ub, maxfes, solverOptions);
		
		allfes(:, iruns, ifitfun) = out.fes;
		allfvals(:, iruns, ifitfun) = out.fmin;
		allout{iruns, ifitfun} = out;
	end
end

% Computational complexity
fprintf('Computing T0.....');
startT0 = tic;
for i = 1 : 1000000
	x = 0.55 + double(i);
	x = x + x;
	x = x ./ 2;
	x = x * x;
	x = sqrt(x);
	x = log(x);
	x = exp(x);
	y = x/x;
end
T0 = toc(startT0) + y - y;
fprintf('Done\n');
fprintf('Computing T1.....');
startT1 = tic;
for i = 1 : 200000
	feval('cec13_f14', rand(D, 1));
end
T1 = toc(startT1);
fprintf('Done\n');
fprintf('Computing T2');
allT2 = zeros(5, 1);
solverOptions.RecordPoint = 0;
for i = 1 : 5
	fprintf('.');
	startT2 = tic;
	feval(solver, 'cec13_f14', lb, ub, 200000, solverOptions);
	allT2(i) = toc(startT2);
end
T2 = mean(allT2);
fprintf('Done\n');
end