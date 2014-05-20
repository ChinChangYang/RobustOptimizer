function [allout, allfvals, allfes] = complete_cec11(...
	solver, ...
	measureOptions, ...
	solverOptions)
% COMPLETE_CEC11 complete four of CEC'11 experiments

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

load('InitialX_CEC11.mat');
defaultMeasureOptions.Runs = 2;
defaultMeasureOptions.FitnessFunctions = ...
	{'cec11_f1', 'cec11_f7', 'cec11_f12', 'cec11_f13'};
measureOptions = setdefoptions(measureOptions, defaultMeasureOptions);
runs	= measureOptions.Runs;
fitfuns = measureOptions.FitnessFunctions;

defaultSolverOptions.dimensionFactor = 5;
defaultSolverOptions.RecordPoint = 21;
defaultSolverOptions.ftarget = -inf;
defaultSolverOptions.TolX = 0;
defaultSolverOptions.TolStagnationIteration = Inf;

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
	
	if isequal(fitfun, 'cec11_f1')
		D		= 6;
		maxfes	= D * 1e4;
		lb		= -6.4 * ones(D, 1);
		ub		= 6.35 * ones(D, 1);	
		solverOptions.NP = 5 * D;
		solverOptions.initial.X = eval(sprintf('XF%dNP%d', ...
			1, ...
			solverOptions.NP));
	elseif isequal(fitfun, 'cec11_f7')
		D		= 20;
		maxfes	= D * 1e4;
		lb		= zeros(D, 1);
		ub		= 2 * pi * ones(D, 1);	
		solverOptions.NP = 5 * D;	
		solverOptions.initial.X = eval(sprintf('XF%dNP%d', ...
			2, ...
			solverOptions.NP));
	elseif isequal(fitfun, 'cec11_f12')
		D		= 26;
		maxfes	= D * 1e4;
		[lb, ub] = getlimit_messenger;
		solverOptions.NP = 5 * D;	
		solverOptions.initial.X = eval(sprintf('XF%dNP%d', ...
			3, ...
			solverOptions.NP));
	elseif isequal(fitfun, 'cec11_f13')
		D		= 22;
		maxfes	= D * 1e4;
		[lb, ub] = getlimit_cassini2;
		solverOptions.NP = 5 * D;	
		solverOptions.initial.X = eval(sprintf('XF%dNP%d', ...
			4, ...
			solverOptions.NP));
	else
		error('BUG: Unknown fitfun');
	end
	
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

fprintf('Done\n');
end