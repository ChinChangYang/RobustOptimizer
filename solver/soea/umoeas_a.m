function [xmin, fmin, out] = umoeas_a(fitfun, lb, ub, maxfunevals, options)
% UMOEAS_A United Multi-Operator Evolutionary Algorithms (variant A)
% UMOEAS_A(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% UMOEAS_A(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP1 = 600;
defaultOptions.NP2 = 100;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.usefunevals = inf;
defaultOptions.initial.X1 = [];
defaultOptions.initial.X2 = [];
defaultOptions.initial.fx1 = [];
defaultOptions.initial.fx2 = [];
defaultOptions.initial.counteval = [];
defaultOptions.initial.countiter = [];
defaultOptions.initial.countcon = [];
defaultOptions.nonlcon = [];
defaultOptions.EarlyStop = 'none';

options = setdefoptions(options, defaultOptions);
usefunevals = options.usefunevals;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

if ~isempty(strfind(options.EarlyStop, 'fitness'))
	EarlyStopOnFitness = true;
	AutoEarlyStop = false;
elseif ~isempty(strfind(options.EarlyStop, 'auto'))
	EarlyStopOnFitness = false;
	AutoEarlyStop = true;
else
	EarlyStopOnFitness = false;
	AutoEarlyStop = false;
end

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X1		= options.initial.X1;
	X2		= options.initial.X2;
	fx1		= options.initial.fx1;
	fx2		= options.initial.fx2;
	counteval = options.initial.counteval;
	countiter = options.initial.countiter;
else
	X1		= [];
	X2		= [];
	fx1		= [];
	fx2		= [];
	counteval = [];
	countiter = [];
end

D = numel(lb);

if isempty(X1)
	NP1 = options.NP1;
else
	[~, NP1] = size(X1);
end

if isempty(X2)
	NP2 = options.NP2;
else
	[~, NP2] = size(X2);
end

% Initialize variables
out = initoutput(RecordPoint, D, NP1, maxfunevals);

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

% Initialize population
if isempty(X1)
	X1 = zeros(D, NP1);
	for i = 1 : NP1
		X1(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
end
if isempty(X2)
	X2 = zeros(D, NP2);
	for i = 1 : NP2
		X2(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
end

% Evaluation
if isempty(fx1)
	fx1 = zeros(1, NP1);
	for i = 1 : NP1
		fx1(i) = feval(fitfun, X1(:, i));
		counteval = counteval + 1;
	end
end
if isempty(fx2)
	fx2 = zeros(1, NP2);
	for i = 1 : NP2
		fx2(i) = feval(fitfun, X2(:, i));
		counteval = counteval + 1;
	end
end

% Sort
[fx1, fidx1] = sort(fx1);
X1 = X1(:, fidx1);
[fx2, fidx2] = sort(fx2);
X2 = X2(:, fidx2);

% Initialize variables
bestSolver = 1;
if D == 10
	CS = 50;
else
	CS = 100;
end
success = zeros(CS, 2);
if fx1(1) < fx2(1)
	fmin = fx1(1);
	xmin = X1(:, 1);
else
	fmin = fx2(1);
	xmin = X2(:, 1);
end
mu = floor(0.5 * NP2);
w = log(mu + 0.5) - log(1 : mu)';
w = w / sum(w);
resumeOptions1 = options;
resumeOptions1.NP = options.NP1;
resumeOptions1.EarlyStop = 'none';
resumeOptions1.usefunevals = max(NP1, NP2);
resumeOptions1.RecordPoint = 0;
resumeOptions1.initial.counteval = 0;
resumeOptions2 = options;
resumeOptions2.NP = options.NP2;
resumeOptions2.EarlyStop = 'none';
resumeOptions2.usefunevals = max(NP1, NP2);
resumeOptions2.RecordPoint = 0;
resumeOptions2.initial.counteval = 0;
out1.final.counteval = 0;
out2.final.counteval = 0;

% Display
if isDisplayIter
	displayitermessages(...
		X1, X2, fx1, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, [X1, X2], [fx1, fx2], counteval, countiter);

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP1;
	outofusefunevals = counteval > usefunevals - NP1;
	if ~EarlyStopOnFitness && ~AutoEarlyStop
		if outofmaxfunevals || outofusefunevals
			break;
		end
	elseif AutoEarlyStop
		X = [X1, X2];
		reachftarget = min([fx1, fx2]) <= ftarget;
		TolX = 10 * eps(mean(X(:)));
		solutionconvergence = std(X(:)) <= TolX;
		TolFun = 10 * eps(mean([fx1, fx2]));
		functionvalueconvergence = std([fx1, fx2]) <= TolFun;
		stagnation = countstagnation >= TolStagnationIteration;
		
		if outofmaxfunevals || ...
				reachftarget || ...
				solutionconvergence || ...
				functionvalueconvergence || ...
				stagnation
			break;
		end
	elseif EarlyStopOnFitness
		reachftarget = fmin <= ftarget;
		
		if outofmaxfunevals || ...
				reachftarget
			break;
		end
	end
	
	% Iteration counter
	countiter = countiter + 1;
	
	if countiter == CS + 1
		% Decide which multi-operator algorithm is the best.
		expected = zeros(1, 2);
		
		for i = 1 : 2
			si= ceil(CS / 2) : CS;
			p = fit(si', success(si, i), 'exp1');
			expected(i) = (p.a * exp(p.b * 2 * CS));
		end
		
		[~, bestSolver] = min(expected);
	elseif countiter == 2 * CS
		% Calculate the mean(x) and standard deviation (sigma) vector of
		% the \mu best individuals of the best_EA, and replace the k-th
		% individual with xk = N(x, sigma), where k is the second worst
		% individual.
		% Replace the worst individuals in the worst performing
		% multi-operator algorithms by the best individual found so far.
		if min(fx1) > fmin
			X1(:, end) = xmin;
			fx1(end) = fmin;
		end
		
		if min(fx2) > fmin
			X2(:, end) = xmin;
			fx2(end) = fmin;
		end
		
		if bestSolver == 1
			xmean1 = mean(X1(:, 1 : 2), 2);
			xstd1 = std(X1(:, 1 : 2), [], 2);
			X2(:, end - 1) = xmean1 + xstd1 .* randn(D, 1);
			fx2(:, end - 1) = feval(fitfun, X2(:, end - 1));
			counteval = counteval + 1;
		elseif bestSolver == 2
			xmean2 = mean(X2(:, 1 : 2), 2);
			xstd2 = std(X2(:, 1 : 2), [], 2);
			X1(:, end - 1) = xmean2 + xstd2 .* randn(D, 1);
			fx1(:, end - 1) = feval(fitfun, X1(:, end - 1));
			counteval = counteval + 1;
		end
		
		[fx1, fidx1] = sort(fx1);
		X1 = X1(:, fidx1);
		resumeOptions1.initial.X = X1;
		resumeOptions1.initial.f = fx1;
		
		[fx2, fidx2] = sort(fx2);
		X2 = X2(:, fidx2);
		resumeOptions2.initial.m = X2(:, 1 : mu) * w;
		
		countiter = 1;
	end
	
	if countiter < CS || bestSolver == 1
		% Evolve population 1 using multi-operator DE.
		[~, ~, out1] = ...
			lshade_sps_eig_j(...
			fitfun, ...
			lb, ...
			ub, ...
			maxfunevals - counteval + out1.final.counteval, ...
			resumeOptions1);
		
		X1 = out1.final.X;
		fx1 = out1.final.f;
		counteval = ...
			counteval ...
			+ out1.final.counteval ...
			- resumeOptions1.initial.counteval;
		
		resumeOptions1.initial = out1.final;
	end
	
	if countiter < CS || bestSolver == 2
		% Evolve population 2 using multi-operator ES.
		[~, ~, out2] = ...
			cmaes(...
			fitfun, ...
			lb, ...
			ub, ...
			maxfunevals - counteval + out2.final.counteval, ...
			resumeOptions2);
		
		X2 = out2.final.X;
		fx2 = out2.final.f;
		counteval = ...
			counteval ...
			+ out2.final.counteval ...
			- resumeOptions2.initial.counteval;
		
		resumeOptions2.initial = out2.final;
	end
	
	resumeOptions1.usefunevals = ...
		out1.final.counteval + max(numel(out1.final.f), numel(out2.final.f));
	resumeOptions2.usefunevals = ...
		out2.final.counteval + max(numel(out1.final.f), numel(out2.final.f));
	
	if countiter <= CS
		success(countiter, :) = [fx1(1), fx2(1)];
	end
	
	if fx1(1) < fx2(1)
		fmin = fx1(1);
		xmin = X1(:, 1);
	else
		fmin = fx2(1);
		xmin = X2(:, 1);
	end
	
	% Record
	out = updateoutput(out, [X1, X2], [fx1, fx2], counteval, countiter);
end

final.X1 = X1;
final.X2 = X2;
final.fx1 = fx1;
final.fx2 = fx2;
final.counteval = counteval;
final.countiter = countiter;

out = finishoutput(out, [X1, X2], [fx1, fx2], counteval, countiter, ...
	'final', final);
end
