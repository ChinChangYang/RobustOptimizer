function [xmin, fmin, out] = sadebin(fitfun, lb, ub, maxfunevals, options)
% SADEBIN SaDE
% SADEBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% SADEBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.LP = 40;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.CRMemory = [];
defaultOptions.initial.ns = [];
defaultOptions.initial.nf = [];
defaultOptions.initial.g = [];
defaultOptions.ConstraintHandling = 'Interpolation';

options = setdefoptions(options, defaultOptions);
LP = round(options.LP);
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

if isequal(options.ConstraintHandling, 'Interpolation')
	interpolation = true;
else
	interpolation = false;
end

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	fx = options.initial.f;
	CRMemory = options.initial.CRMemory;
	ns = options.initial.ns;
	nf = options.initial.nf;
	g = options.initial.g;
else
	X = [];
	fx = [];
	CRMemory = [];
	ns = [];
	nf = [];
	g = [];
end

D = numel(lb);
if isempty(X)	
	NP = options.NP;
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'FC');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
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

% Initialize variables
K = 4;

% CRMemory
if isempty(CRMemory)
	CRMemory = cell(1, K);
	for i = 1 : K
		CRMemory{i} = 0.5 + 0.1 * randn;
	end
end

% ns
if isempty(ns)
	ns = ones(K, LP);
end

% nf
if isempty(nf)
	nf = ones(K, LP);
end

% g
if isempty(g)
	g = 2;
end

V = X;
U = X;
selStrategy = zeros(1, NP);
CR = zeros(K, NP);
CRm = 0.5 * ones(K, 1);
p = 1/K * ones(K, 1);
fu = zeros(1, NP);
FC = zeros(1, NP);		% Consecutive Failure Counter
rt = zeros(1, NP);
r1 = zeros(1, NP);
r2 = zeros(1, NP);
r3 = zeros(1, NP);
r4 = zeros(1, NP);
r5 = zeros(1, NP);

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter, ...
	'FC', FC);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(fx) <= ftarget;
	stagnation = countStagnation >= TolStagnationIteration;	
	if outofmaxfunevals || reachftarget || stagnation
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
	
	for i = 1 : NP
		rt(i) = i;
		
		% Generate r1
		r1(i) = floor(1 + NP * rand);
		while rt(i) == r1(i)
			r1(i) = floor(1 + NP * rand);
		end
		
		% Generate r2
		r2(i) = floor(1 + NP * rand);
		while rt(i) == r2(i) || r1(i) == r2(i)
			r2(i) = floor(1 + NP * rand);
		end
		
		% Generate r3
		r3(i) = floor(1 + NP * rand);
		while rt(i) == r3(i) || r1(i) == r3(i) || r2(i) == r3(i)
			r3(i) = floor(1 + NP * rand);
		end
		
		% Generate r4
		r4(i) = floor(1 + NP * rand);
		while rt(i) == r4(i) || r1(i) == r4(i) || r2(i) == r4(i) ...
				|| r3(i) == r4(i)
			r4(i) = floor(1 + NP * rand);
		end
		
		% Generate r5
		r5(i) = floor(1 + NP * rand);
		while rt(i) == r5(i) || r1(i) == r5(i) || r2(i) == r5(i) ...
				|| r3(i) == r5(i) || r4(i) == r5(i)
			r5(i) = floor(1 + NP * rand);
		end
	end
	
	% Mutation
	[~, besti] = min(fx);
	for i = 1 : NP
		if selStrategy(rt(i)) == 1
			% DE/rand/1			
			V(:, i) = X(:, r1(i)) + F(rt(i)) * (X(:, r2(i)) - X(:, r3(i)));
		elseif selStrategy(rt(i)) == 2
			% DE/current-to-best/2
			V(:, i) = X(:, rt(i)) + F(rt(i)) * (X(:, besti) - X(:, rt(i))) + ...
				F(i) * (X(:, r1(i)) - X(:, r2(i))) + ...
				F(i) * (X(:, r3(i)) - X(:, r4(i)));
		elseif selStrategy(rt(i)) == 3
			% DE/rand/2
			V(:, i) = X(:, r1(i)) + F(rt(i)) * (X(:, r2(i)) - X(:, r3(i))) + ...
				F(i) * (X(:, r4(i)) - X(:, r5(i)));
		elseif selStrategy(rt(i)) == 4
			% DE/current-to-rand/1
			V(:, i) = X(:, rt(i)) + rand * (X(:, r1(i)) - X(:, rt(i))) + ...
				F(rt(i)) * (X(:, r2(i)) - X(:, r3(i)));
		else
			fprintf('BUG\n');
		end
	end
	
	for i = 1 : NP
		if selStrategy(rt(i)) == 1
			% DE/rand/1/bin
			jrand = floor(1 + D * rand);
			for j = 1 : D
				if rand < CR(1, rt(i)) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, rt(i));
				end
			end
		elseif selStrategy(rt(i)) == 2
			% DE/rand-to-best/1/bin
			jrand = floor(1 + D * rand);
			for j = 1 : D
				if rand < CR(2, rt(i)) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, rt(i));
				end
			end
		elseif selStrategy(rt(i)) == 3
			% DE/rand/2/bin
			jrand = floor(1 + D * rand);
			for j = 1 : D
				if rand < CR(3, rt(i)) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, rt(i));
				end
			end
		elseif selStrategy(rt(i)) == 4
			% DE/current-to-rand/1
			U(:, i) = V(:, i);
		else
			fprintf('BUG\n');
		end
	end
		
	if interpolation
		% Correction for outside of boundaries
		for i = 1 : NP
			for j = 1 : D
				if U(j, i) < lb(j)
					U(j, i) = 0.5 * (lb(j) + X(j, rt(i)));
				elseif U(j, i) > ub(j)
					U(j, i) = 0.5 * (ub(j) + X(j, rt(i)));
				end
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			X, U, fx, countiter, XX, YY, ZZ);
	end
	
	% Evaluation
	for i = 1 : NP
		fu(i) = feval(fitfun, U(:, i));
		counteval = counteval + 1;
	end
	
	% Selection
	FailedIteration = true;
	for i = 1 : NP
		k = selStrategy(i);		
		if fu(i) < fx(i)
			X(:, i)				= U(:, i);
			fx(i)				= fu(i);
			FailedIteration		= false;
			FC(i)				= 0;
			CRMemory{k}(end + 1) = CR(k, rt(i));
			if g <= LP
				ns(k, g) = ns(k, g) + 1;
			else
				ns(k, end) = ns(k, end) + 1;
			end
		else
			FC(i) = FC(i) + 1;
			if g <= LP
				nf(k, g) = nf(k, g) + 1;
			else
				nf(k, end) = nf(k, end) + 1;
			end
		end
	end
	
	% Record
	out = updateoutput(out, X, fx, counteval, countiter, ...
		'FC', FC);
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end	
end

[fmin, minindex] = min(fx);
xmin = X(:, minindex);

out = finishoutput(out, X, fx, counteval, countiter, ...
	'FC', zeros(NP, 1));
end