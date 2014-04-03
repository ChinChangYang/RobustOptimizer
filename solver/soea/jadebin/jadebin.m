function [xmin, fmin, out] = jadebin(fitfun, lb, ub, maxfunevals, options)
% JADEBIN JADE algorithm
% JADEBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% JADEBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.F = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.delta_CR = 0.1;
defaultOptions.delta_F = 0.1;
defaultOptions.p = 0.05;
defaultOptions.w = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.mu_CR = [];
defaultOptions.initial.mu_F = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
delta_CR = options.delta_CR;
delta_F = options.delta_F;
p = options.p;
w = options.w;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
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
	NP = max(1, ceil(dimensionFactor * D));
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'MF', 'MCR');

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
	f = zeros(1, NP);
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
end

% Sort
[f, fidx] = sort(f);
X = X(:, fidx);

% mu_F
if isempty(mu_F)
	mu_F = options.F;
end

% mu_CR
if isempty(mu_CR)
	mu_CR = options.CR;
end

% Initialize variables
V = X;
U = X;
pbest_size = p * NP;
A_size = 0;

% Display
if isDisplayIter
	displayitermessages(...
		X, U, f, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, f, counteval, ...
	'MF', mu_F, 'MCR', mu_CR);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(f) <= ftarget;
	stagnation = countStagnation >= TolStagnationIteration;
	
	% Convergence conditions	
	if outofmaxfunevals || reachftarget || stagnation
		break;
	end
	
	% Scaling factor and crossover rate
	nS = 0;
	S_F = zeros(1, NP);
	S_CR = zeros(1, NP);
	CR = mu_CR + delta_CR * randn(1, NP);
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	F = cauchyrnd(mu_F, delta_F, NP, 1);
	F(F > 1) = 1;
	
	while any(F <= 0)		
		F(F <= 0) = cauchyrnd(mu_F, delta_F, sum(F <= 0), 1);
		F(F > 1) = 1;
	end
	
	XA = [X, A];
	
	% Mutation
	for i = 1 : NP				
		% Generate pbest_idx
		pbest_idx = max(1, ceil(rand * pbest_size));
		
		% Generate r1
		r1 = floor(1 + NP * rand);
		while i == r1
			r1 = floor(1 + NP * rand);
		end
		
		% Generate r2
		r2 = floor(1 + (NP + A_size) * rand);
		while i == r1 || r1 == r2
			r2 = floor(1 + (NP + A_size) * rand);
		end
		
		V(:, i) = X(:, i) + F(i) .* ...
			(X(:, pbest_idx) - X(:, i) + X(:, r1) - XA(:, r2));
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
	
	% Correction for outside of boundaries
	for i = 1 : NP
		for j = 1 : D
			if U(j, i) < lb(j)
				U(j, i) = 0.5 * (lb(j) + X(j, i));
			elseif U(j, i) > ub(j)
				U(j, i) = 0.5 * (ub(j) + X(j, i));
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			X, U, f, countiter, XX, YY, ZZ);
	end
	
	% Selection
	FailedIteration = true;
	for i = 1 : NP
		fui = feval(fitfun, U(:, i));
		counteval = counteval + 1;
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
			S_CR(nS + 1) = CR(i);
			S_F(nS + 1) = F(i);
			nS = nS + 1;
			FailedIteration = false;
			
			if A_size < NP
				A_size = A_size + 1;
				A(:, A_size) = X(:, i);
			else
				ri = floor(1 + NP * rand);
				A(:, ri) = X(:, i);
			end
		end
	end
	
	% Update CR and F
	if nS > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : nS));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : nS).^2) / sum(S_F(1 : nS));
	end
	
	% Sort	
	[f, fidx] = sort(f);
	X = X(:, fidx);
	
	% Record
	out = updateoutput(out, X, f, counteval, ...
		'MF', mu_F, 'MCR', mu_CR);
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end	
end

fmin = f(1);
xmin = X(:, 1);

if fmin < out.bestever.fmin
	out.bestever.fmin = fmin;
	out.bestever.xmin = xmin;
end

final.A = A;
final.mu_F = mu_F;
final.mu_CR = mu_CR;

out = finishoutput(out, X, f, counteval, ...
	'final', final, ...
	'MF', mu_F, ...
	'MCR', mu_CR);
end
