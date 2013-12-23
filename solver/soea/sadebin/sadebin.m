function [xmin, fmin, out] = sadebin(fitfun, lb, ub, maxfunevals, options)
% SADEBIN SaDE
% SADEBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% SADEBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.maxfunevalsFactor = 0;
defaultOptions.LP = 40;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.CRMemory = [];
defaultOptions.initial.ns = [];
defaultOptions.initial.nf = [];
defaultOptions.initial.g = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
LP = round(options.LP);
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
D = numel(lb);
X = options.initial.X;
f = options.initial.f;
CRMemory = options.initial.CRMemory;
ns = options.initial.ns;
nf = options.initial.nf;
g = options.initial.g;

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
K = 4;

if isempty(CRMemory)
	CRMemory = cell(1, K);
	for i = 1 : K
		CRMemory{i} = 0.5 + 0.1 * randn;
	end
end

if isempty(ns)
	ns = ones(K, LP);
end

if isempty(nf)
	nf = ones(K, LP);
end

if isempty(g)
	g = 2;
end

V = X;
U = X;

selStrategy = zeros(1, NP);
CR = zeros(K, NP);
CRm = 0.5 * ones(K, 1);
p = 1/K * ones(K, 1);

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
	
	if countiter > LP
		S = sum(ns, 2) ./ (sum(ns, 2) + sum(nf, 2)) + eps;
		p = S ./ sum(S);
		ns(:, 1:end-1) = ns(:, 2:end);
		nf(:, 1:end-1) = nf(:, 2:end);
	end
	
	cp = zeros(K, 1);
	cp(1) = p(1);
	for i = 2 : K
		cp(i) = cp(i-1) + p(i);
	end
	
	for i = 1 : NP
		s = rand;
		selStrategy(i) = find(cp >= s, 1, 'first');
	end
	
	F = 0.5 + 0.3 * randn(1, NP);
	
	if g >= LP
		for k = 1 : K
			CRm(k) = median(CRMemory{k});
		end
	end
	
	for k = 1 : K
		CR(k, :) = CRm(k) + 0.1 * randn(1, NP);
		checkIndex = CR(k, :) < 0 | CR(k, :) > 1;
		while any(checkIndex > 0)
			CR(k, checkIndex) = CRm(k) + 0.1 * randn(1, sum(checkIndex));
			checkIndex = CR(k, :) < 0 | CR(k, :) > 1;
		end
	end
	
	% Mutation
	[~, besti] = min(f);
	for i = 1 : NP
		if selStrategy(i) == 1
			% DE/rand/1
			r1 = floor(1 + NP * rand);
			r2 = floor(1 + NP * rand);
			r3 = r2;
			
			while r2 == r3
				r3 = floor(1 + NP * rand);
			end
			
			pF = 1 + 0.01 * randn;
			V(:, i) = X(:, r1) + F(i) * pF * (X(:, r2) - X(:, r3));
		elseif selStrategy(i) == 2
			% DE/current-to-best/2
			r2 = floor(1 + NP * rand);
			r3 = r2;
			r4 = floor(1 + NP * rand);
			r5 = r4;
			
			while r2 == r3
				r3 = floor(1 + NP * rand);
			end
			
			while r4 == r5
				r5 = floor(1 + NP * rand);
			end
			
			pF1 = 1 + 0.03 * randn;
			pF2 = 1 + 0.03 * randn;
			pF3 = 1 + 0.03 * randn;
			V(:, i) = X(:, i) + F(i) * pF1 * (X(:, besti) - X(:, i)) + ...
				F(i) * pF2 * (X(:, r2) - X(:, r3)) + ...
				F(i) * pF3 * (X(:, r4) - X(:, r5));
		elseif selStrategy(i) == 3
			% DE/rand/2
			r1 = floor(1 + NP * rand);
			r2 = floor(1 + NP * rand);
			r3 = r2;
			r4 = floor(1 + NP * rand);
			r5 = r4;
			
			while r2 == r3
				r3 = floor(1 + NP * rand);
			end
			
			while r4 == r5
				r5 = floor(1 + NP * rand);
			end
			
			pF1 = 1 + 0.03 * randn;
			pF2 = 1 + 0.03 * randn;
			V(:, i) = X(:, r1) + F(i) * pF1 * (X(:, r2) - X(:, r3)) + ...
				F(i) * pF2 * (X(:, r4) - X(:, r5));
		elseif selStrategy(i) == 4
			% DE/current-to-rand/1
			r1 = floor(1 + NP * rand);
			r2 = floor(1 + NP * rand);
			r3 = r2;
			
			while r2 == r3
				r3 = floor(1 + NP * rand);
			end
			
			pF1 = 1 + 0.03 * randn;
			pF2 = 1 + 0.03 * randn;
			V(:, i) = X(:, i) + rand * pF1 * (X(:, r1) - X(:, i)) + ...
				F(i) * pF2 * (X(:, r2) - X(:, r3));
		else
			fprintf('BUG\n');
		end
	end
	
	for i = 1 : NP
		if selStrategy(i) == 1
			% DE/rand/1/bin
			jrand = floor(1 + D * rand);
			for j = 1 : D
				if rand < CR(1, i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, i);
				end
			end
		elseif selStrategy(i) == 2
			% DE/rand-to-best/1/bin
			jrand = floor(1 + D * rand);
			for j = 1 : D
				if rand < CR(2, i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, i);
				end
			end
		elseif selStrategy(i) == 3
			% DE/rand/2/bin
			jrand = floor(1 + D * rand);
			for j = 1 : D
				if rand < CR(3, i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, i);
				end
			end
		elseif selStrategy(i) == 4
			% DE/current-to-rand/1
			U(:, i) = V(:, i);
		else
			fprintf('BUG\n');
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
	for i = 1 : NP
		k = selStrategy(i);
		fui = feval(fitfun, U(:, i));
		counteval = counteval + 1;
		
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
			CRMemory{k}(end + 1) = CR(k, i);
			if g <= LP
				ns(k, g) = ns(k, g) + 1;
			else
				ns(k, end) = ns(k, end) + 1;
			end
		else
			if g <= LP
				nf(k, g) = nf(k, g) + 1;
			else
				nf(k, end) = nf(k, end) + 1;
			end
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

% Final state of this algorithm
final.CRMemory = CRMemory;
final.ns = ns;
final.nf = nf;
final.g = g;
out = finishoutput(out, X, f, counteval, 'final', final);
end