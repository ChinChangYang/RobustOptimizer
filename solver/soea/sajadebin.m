function [xmin, fmin, out] = sajadebin(fitfun, lb, ub, maxfunevals, options)
% JADEBIN Self-adaptive JADE algorithm
% JADEBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% JADEBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.delta_CR = 0.1;
defaultOptions.delta_F = 0.1;
defaultOptions.w = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = 0;
defaultOptions.TolX = 0;
defaultOptions.TolStagnationIteration = 30;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.mu_CR = [];
defaultOptions.initial.mu_F = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
delta_CR = options.delta_CR;
delta_F = options.delta_F;
w = options.w;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
TolStagnationIteration = options.TolStagnationIteration;

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	f = options.initial.f;
	A = options.initial.A;
	mu_CR = options.initial.mu_CR;
	mu_F = options.initial.mu_F;
else
	X = [];
	f = [];
	A = [];
	mu_CR = [];
	mu_F = [];
end

D = numel(lb);
if isempty(X)
	NP = max(4, ceil(dimensionFactor * D));
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'mu_F', 'mu_CR');

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

% Initialize archive
if isempty(A)
	A = X;
end

% Evaluation
if isempty(f)
	f = nan(1, NP);
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
end

[~, pfidx] = sort(f);
f(1 : NP) = f(pfidx);
X(:, 1 : NP) = X(:, pfidx);

% mu_F
if isempty(mu_F)
	mu_F = options.F;
end

% mu_CR
if isempty(mu_CR)
	mu_CR = options.CR;
end

% Initialize variables
r1 = zeros(1, NP);
r2 = zeros(1, NP);
V = X(:, 1 : NP);
U = X(:, 1 : NP);
P1 = ones(NP, NP);
P2 = ones(NP, 2 * NP);

% Display
if isDisplayIter
	displayitermessages(...
		X(:, 1 : NP), U(:, 1 : NP), f(1 : NP), countiter, XX, YY, ZZ, ...
		'mu_F', mu_F, 'mu_CR', mu_CR);
end

% Record
out = updateoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
	'mu_F', mu_F, 'mu_CR', mu_CR);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(f(1 : NP)) <= ftarget;
	fitnessconvergence = isConverged(f(1 : NP), TolFun);
	solutionconvergence = isConverged(X(:, 1 : NP), TolX);
	stagnation = countStagnation >= TolStagnationIteration;
	
	% Convergence conditions	
	if outofmaxfunevals || reachftarget || fitnessconvergence || ...
			solutionconvergence || stagnation
		break;
	end
	
	% Scaling factor and crossover rate
	S_F = zeros(1, NP);
	S_CR = zeros(1, NP);
	CR = mu_CR + delta_CR * randn(1, NP);
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	F = cauchyrnd(mu_F, delta_F, NP, 1);
	F(F > 1) = 1;
	
	for retry = 1 : 3
		if all(F > 0)
			break;
		end
		
		F(F <= 0) = cauchyrnd(mu_F, delta_F, sum(F <= 0), 1);
		F(F > 1) = 1;
	end
	
	F(F <= 0) = 0.01 * mu_F * (1 + rand);
	
	A_Counter = 0;
	XA = [X(:, 1 : NP), A];
	
	% Probability table
	for i = 1 : NP
		P1(i, i) = 0;
		P2(i, i) = 0;
		P1(i, :) = P1(i, :) / sum(P1(i, :));
		P2(i, :) = P2(i, :) / sum(P2(i, :));
	end
	
	% Mutation
	[~, besti] = min(f);
	for i = 1 : NP
		ptable1 = cumsum([0, P1(i, :)]);
		ptable2 = cumsum([0, P2(i, :)]);
		for checkRetry = 1 : 3			
			% Generate r1
			for retry = 1 : 3
				r1(i) = sum(rand >= ptable1);
				if i ~= r1(i)
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : 3
				r2(i) = sum(rand >= ptable2);
				if ~(all(X(:, i) == XA(:, r2(i))) || ...
						all(X(:, r1(i)) == XA(:, r2(i))))
					break;
				end
			end
			
			V(:, i) = X(:, i) + F(i) .* ...
				(X(:, besti) - X(:, i) + X(:, r1(i)) - XA(:, r2(i)));
			
			% Check boundary
			if all(V(:, i) > lb) && all(V(:, i) < ub)
				break;
			end
		end
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		for j = 1 : D
			if rand < CR(i) || j == jrand
				U(j, i) = V(j, i);
			else
				U(j, i) = X(j, i);
			end
		end
	end
	
	% Bounds reflection
	for i = 1 : NP
		for j = 1 : D
			for k = 1 : 3
				if U(j, i) < lb(j)
					U(j, i) = 2 * lb(j) - U(j, i);
				end
				
				if U(j, i) > ub(j)
					U(j, i) = 2 * ub(j) - U(j, i);
				end
				
				if U(j, i) >= lb(j) && U(j, i) <= ub(j)
					break;
				end
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			X(:, 1 : NP), U(:, 1 : NP), f(1 : NP), countiter, XX, YY, ZZ, ...
			'mu_F', mu_F, 'mu_CR', mu_CR);
	end
	
	% Selection
	FailedIteration = true;
	for i = 1 : NP
		fui = feval(fitfun, U(:, i));
		counteval = counteval + 1;
		
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
			A(:, NP + A_Counter + 1) = U(:, i);
			S_CR(A_Counter + 1) = CR(i);
			S_F(A_Counter + 1) = F(i);
			A_Counter = A_Counter + 1;
			FailedIteration = false;
			P1(i, r1(i)) = P1(i, r1(i)) * 2;
			P2(i, r2(i)) = P2(i, r2(i)) * 2;
		else
			P1(i, r1(i)) = P1(i, r1(i)) / 2;
			P2(i, r2(i)) = P2(i, r2(i)) / 2;
		end
	end
	
	% Update archive
	rand_idx = randperm(NP + A_Counter);
	A(:, 1 : NP) = A(:, rand_idx(1 : NP));
	
	% Update CR and F
	if A_Counter > 0
		mu_CR = (1 - w) * mu_CR + w * mean(S_CR(1 : A_Counter));
		mu_F = (1 - w) * mu_F + w * sum(S_F(1 : A_Counter).^2) / sum(S_F(1 : A_Counter));
	end
	
	% Record
	out = updateoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
		'mu_F', mu_F, 'mu_CR', mu_CR);
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end	
end

[fmin, besti] = min(f);
xmin = X(:, besti);

if fmin < out.bestever.fmin
	out.bestever.fmin = fmin;
	out.bestever.xmin = xmin;
end

final.A = A;
final.mu_F = mu_F;
final.mu_CR = mu_CR;

out = finishoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
	'final', final, ...
	'mu_F', mu_F, 'mu_CR', mu_CR);
end
