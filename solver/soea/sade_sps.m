function [xmin, fmin, out] = sade_sps(fitfun, lb, ub, maxfunevals, options)
% SADE_SPS SaDE with SPS Framework
% SADE_SPS(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% SADE_SPS(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.LP = 40;
defaultOptions.Q = 70;
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

options = setdefoptions(options, defaultOptions);
LP = round(options.LP);
Q = options.Q;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

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
SP = X;
fSP = fx;
iSP = 1;

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
	
	% Mutation
	[~, ibestX] = min(fx);
	[~, ibestSP] = min(fSP);
	for i = 1 : NP		
		% Generate r1
		r1 = floor(1 + NP * rand);
		while i == r1
			r1 = floor(1 + NP * rand);
		end
		
		% Generate r2
		r2 = floor(1 + NP * rand);
		while i == r2 || r1 == r2
			r2 = floor(1 + NP * rand);
		end
		
		% Generate r3
		r3 = floor(1 + NP * rand);
		while i == r3 || r1 == r3|| r2 == r3
			r3 = floor(1 + NP * rand);
		end
		
		% Generate r4
		r4 = floor(1 + NP * rand);
		while i == r4 || r1 == r4 || r2 == r4 || r3 == r4
			r4 = floor(1 + NP * rand);
		end
		
		% Generate r5
		r5 = floor(1 + NP * rand);
		while i == r5 || r1 == r5 || r2 == r5 || r3 == r5 || r4 == r5
			r5 = floor(1 + NP * rand);
		end
		
		if FC(i) <= Q
			if selStrategy(i) == 1
				% DE/rand/1
				V(:, i) = X(:, r1) + F(i) * (X(:, r2) - X(:, r3));
			elseif selStrategy(i) == 2
				% DE/current-to-best/2
				V(:, i) = X(:, i) + F(i) * (X(:, ibestX) - X(:, i)) + ...
					F(i) * (X(:, r1) - X(:, r2)) + ...
					F(i) * (X(:, r3) - X(:, r4));
			elseif selStrategy(i) == 3
				% DE/rand/2
				V(:, i) = X(:, r1) + F(i) * (X(:, r2) - X(:, r3)) + ...
					F(i) * (X(:, r4) - X(:, r5));
			elseif selStrategy(i) == 4
				% DE/current-to-rand/1
				V(:, i) = X(:, i) + rand * (X(:, r1) - X(:, i)) + ...
					F(i) * (X(:, r2) - X(:, r3));
			else
				fprintf('ASSERTION: SHOULD NEVER REACHED HERE!\n');
			end
		else
			if selStrategy(i) == 1
				% DE/rand/1
				V(:, i) = SP(:, r1) + F(i) * (SP(:, r2) - SP(:, r3));
			elseif selStrategy(i) == 2
				% DE/current-to-best/2
				V(:, i) = SP(:, i) + F(i) * (SP(:, ibestSP) - SP(:, i)) + ...
					F(i) * (SP(:, r1) - SP(:, r2)) + ...
					F(i) * (SP(:, r3) - SP(:, r4));
			elseif selStrategy(i) == 3
				% DE/rand/2
				V(:, i) = SP(:, r1) + F(i) * (SP(:, r2) - SP(:, r3)) + ...
					F(i) * (SP(:, r4) - SP(:, r5));
			elseif selStrategy(i) == 4
				% DE/current-to-rand/1
				V(:, i) = SP(:, i) + rand * (SP(:, r1) - SP(:, i)) + ...
					F(i) * (SP(:, r2) - SP(:, r3));
			else
				fprintf('ASSERTION: SHOULD NEVER REACHED HERE!\n');
			end
		end
	end
	
	for i = 1 : NP
		jrand = floor(1 + D * rand);
		if FC(i) <= Q
			if selStrategy(i) == 1
				% DE/rand/1/bin
				for j = 1 : D
					if rand < CR(1, i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = X(j, i);
					end
				end
			elseif selStrategy(i) == 2
				% DE/rand-to-best/1/bin
				for j = 1 : D
					if rand < CR(2, i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = X(j, i);
					end
				end
			elseif selStrategy(i) == 3
				% DE/rand/2/bin
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
				fprintf('ASSERTION: SHOULD NEVER REACHED HERE!\n');
			end
		else			
			if selStrategy(i) == 1
				% DE/rand/1/bin
				for j = 1 : D
					if rand < CR(1, i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = SP(j, i);
					end
				end
			elseif selStrategy(i) == 2
				% DE/rand-to-best/1/bin
				for j = 1 : D
					if rand < CR(2, i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = SP(j, i);
					end
				end
			elseif selStrategy(i) == 3
				% DE/rand/2/bin
				for j = 1 : D
					if rand < CR(3, i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = SP(j, i);
					end
				end
			elseif selStrategy(i) == 4
				% DE/current-to-rand/1
				U(:, i) = V(:, i);
			else
				fprintf('ASSERTION: SHOULD NEVER REACHED HERE!\n');
			end
		end
	end
	
	% Correction for outside of boundaries
	for i = 1 : NP
		if FC(i) <= Q
			for j = 1 : D
				if U(j, i) < lb(j)
					U(j, i) = 0.5 * (lb(j) + X(j, i));
				elseif U(j, i) > ub(j)
					U(j, i) = 0.5 * (ub(j) + X(j, i));
				end
			end
		else
			for j = 1 : D
				if U(j, i) < lb(j)
					U(j, i) = 0.5 * (lb(j) + SP(j, i));
				elseif U(j, i) > ub(j)
					U(j, i) = 0.5 * (ub(j) + SP(j, i));
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
			SP(:, iSP)	= U(:, i);
			fSP(iSP)	= fu(i);
			iSP			= mod(iSP, NP) + 1;
			FailedIteration		= false;
			FC(i)				= 0;
			CRMemory{k}(end + 1) = CR(k, i);
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