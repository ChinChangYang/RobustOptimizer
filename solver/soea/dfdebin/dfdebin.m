function [xmin, fmin, out] = dfdebin(fitfun, lb, ub, maxfunevals, options)
% DFDEBIN Drift-free Differential Evolution
% DFDEBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DFDEBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 2 * numel(lb);
defaultOptions.maxfunevalsFactor =  0;
defaultOptions.F = 1;
defaultOptions.pmu = 0.5;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
F = options.F;
pmu = options.pmu;
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
	for i = 1 : NP
		r1 = floor(1 + NP * rand);
		r2 = floor(1 + NP * rand);
		
		while r1 == i
			r1 = floor(1 + NP * rand);
		end
		
		while r2 == i || r2 == r1
			r2 = floor(1 + NP * rand);
		end
		
		if rand <= pmu
			U(:, i) = X(:, i) + F * (X(:, r1) - X(:, r2));
		else
			sum2 = 0;
			while sum2 == 0
				r3 = floor(1 + NP * rand);
				r4 = floor(1 + NP * rand);
				
				while r3 == i
					r3 = floor(1 + NP * rand);
				end
				
				while r4 == i || r4 == r3
					r4 = floor(1 + NP * rand);
				end
				
				sum1 = 0;
				sum2 = 0;
				
				for j = 1 : D
					d1 = X(j, r1) - X(j, r2);
					d2 = X(j, r3) - X(j, r4) - 2 * X(j, i);
					sum1 = sum1 + d1 * d2;
					sum2 = sum2 + d2 * d2;
				end
			end
			
			K = sqrt(D) * sum1 / sum2;
			
			U(:, i) = X(:, i) + K * (X(:, r3) + X(:, r4) - 2 * X(:, i));
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