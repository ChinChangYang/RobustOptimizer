function [xmin, fmin, out] = decurrenttobest1bin(fitfun, lb, ub, maxfunevals, options)
% DECURRENTTOBEST1BIN Classical DE/current-to-best/1/bin
% DECURRENTTOBEST1BIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DECURRENTTOBEST1BIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.CR = 0.5;
defaultOptions.F = 0.7;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
CR = options.CR;
F = options.F;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
D = numel(lb);
X = options.initial.X;
f = options.initial.f;

if isempty(X)
	NP = ceil(dimensionFactor * D);
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D, NP, maxfunevals);

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = preparecontourdata(D, lb, ub, fitfun);
end

% Initialize population
if isempty(X)
	if NP < 1e1
		LHS = lhsdesign(NP, D, 'iteration', 10)';
	elseif NP < 1e2
		LHS = lhsdesign(NP, D, 'iteration', 2)';
	else
		LHS = rand(D, NP);
	end
	
	X = zeros(D, NP);
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* LHS(:, i);
	end
end

% Evaluation
if isempty(f)
	f = zeros(1, NP);
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
end

% Initialize variables
V = X;
U = X;

% Display
if isDisplayIter
	displayitermessages(X, U, f, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, f, counteval);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions	
	stdf = std(f);
	stdX = std(X, 0, 2);
	meanX = mean(X, 2);
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(f) <= ftarget;
	fitnessconvergence = stdf <= mean(abs(f)) * 100 * eps || stdf <= realmin || stdf <= TolFun;
	solutionconvergence = all(stdX <= meanX * 100 * eps) || all(stdX <= 100 * realmin) || ...
			all(stdX <= TolX);
	
	% Convergence conditions	
	if outofmaxfunevals || reachftarget || fitnessconvergence || solutionconvergence
		break;
	end
	
	% Mutation
	[~, gbest] = min(f);
	for i = 1 : NP
		r1 = floor(1 + NP * rand);
		r2 = floor(1 + NP * rand);
		
		while r1 == r2
			r2 = floor(1 + NP * rand);
		end
		
		V(:, i) = X(:, i) + (F + 0.01 * randn) * (X(:, gbest) - X(:, i)) + ...
			(F + 0.01 * randn) * (X(:, r1) - X(:, r2));
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		for j = 1 : D
			if rand < CR || j == jrand
				U(j, i) = V(j, i);
			else
				U(j, i) = X(j, i);
			end
		end
	end
	
	% Repair
	for i = 1 : NP
		for j = 1 : D
			if U(j, i) < lb(j)
				U(j, i) = X(j, i) + rand * (lb(j) - X(j, i));
			elseif U(j, i) > ub(j)
				U(j, i) = X(j, i) + rand * (ub(j) - X(j, i));
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(X, U, f, countiter, XX, YY, ZZ);
	end
	
	% Selection
	counteval = counteval + NP;
	for i = 1 : NP
		fui = feval(fitfun, U(:, i));
		
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
		end
	end
	
	% Record
	out = updateoutput(out, X, f, counteval);
	
	% Iteration counter
	countiter = countiter + 1;
end

[fmin, fminidx] = min(f);
xmin = X(:, fminidx);

if fmin < out.bestever.fmin
	out.bestever.fmin = fmin;
	out.bestever.xmin = xmin;
end

out = finishoutput(out, X, f, counteval);
end