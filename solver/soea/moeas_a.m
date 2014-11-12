function [xmin, fmin, out] = moeas_a(fitfun, lb, ub, maxfunevals, options)
% MOEAS_A Multi-Operator Evolutionary Algorithms (variant A)
% MOEAS_A(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% MOEAS_A(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP1 = 950;
defaultOptions.NP2 = 150;
defaultOptions.TrialStageFactor = 0.2;
defaultOptions.TQ = 5;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.usefunevals = inf;
defaultOptions.initial.X = [];
defaultOptions.initial.fx = [];
defaultOptions.initial.counteval = [];
defaultOptions.initial.countiter = [];
defaultOptions.initial.countcon = [];
defaultOptions.initial.solver = [];
defaultOptions.nonlcon = [];
defaultOptions.EarlyStop = 'none';

options = setdefoptions(options, defaultOptions);
trialStageFactor = options.TrialStageFactor;
TQ = options.TQ;
isDisplayIter = strcmp(options.Display, 'iter');

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X			= options.initial.X;
	fx			= options.initial.fx;
	counteval	= options.initial.counteval;
	countiter	= options.initial.countiter;
	solver		= options.initial.solver;
else
	X			= [];
	fx			= [];
	counteval	= [];
	countiter	= [];
	solver		= [];
end

D = numel(lb);

if isempty(X)
	NP = options.NP1;
else
	[~, NP] = size(X);
end

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
end

% counteval
if isempty(counteval)
	counteval = 0;
end

% countiter
if isempty(countiter)
	countiter = 1;
end

% solver
if isempty(solver)
	solver = 'lshade_sps';
end

% Initialize population
if isempty(X)
	X = zeros(D, NP);
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
end

% Evaluation
if isempty(fx)
	fx = zeros(1, NP);
	for i = 1 : NP
		fx(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
end

% Sort
[fx, fidx] = sort(fx);
X = X(:, fidx);

% Initialize variables
trialStage = trialStageFactor * maxfunevals;
resumeOptions = options;
resumeOptions.NP = options.NP1;
resumeOptions.EarlyStop = 'fitness';
resumeOptions.usefunevals = trialStage;
resumeOptions.initial.X = X;
resumeOptions.initial.fx = fx;
resumeOptions.initial.counteval = 0;

% Display
if isDisplayIter
	displayitermessages(...
		X, X, fx, countiter, XX, YY, ZZ);
end

[~, ~, innerOut] = ...
	feval(solver, ...
	fitfun, ...
	lb, ...
	ub, ...
	maxfunevals, ...
	resumeOptions);

resumeOptions.initial = innerOut.final;

if innerOut.muFC(end) <= TQ	
	solver = 'cmaes';
	resumeOptions = options;
	resumeOptions.NP = options.NP2;
	resumeOptions.EarlyStop = 'fitness';
	resumeOptions.ConstraintHandling = 'OnBound';
	resumeOptions.initial.counteval = innerOut.final.counteval;
end

resumeOptions.usefunevals = maxfunevals;
[~, ~, innerOut] = ...
	feval(solver, ...
	fitfun, ...
	lb, ...
	ub, ...
	maxfunevals, ...
	resumeOptions);

X = innerOut.final.X;
fx = innerOut.final.f;
fmin = fx(1);
xmin = X(:, 1);
out = innerOut;
end
