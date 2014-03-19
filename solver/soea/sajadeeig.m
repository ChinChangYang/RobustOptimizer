function [xmin, fmin, out] = sajadeeig(fitfun, lb, ub, maxfunevals, options)
% SAJADEEIG Self-adaptive JADE algorithm with EIG crossover
% SAJADEEIG(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% SAJADEEIG(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.F = 0.7;
defaultOptions.H = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.R = 0.5;
defaultOptions.delta_CR = 0.1;
defaultOptions.delta_F = 0.1;
defaultOptions.delta_H = 0.1;
defaultOptions.w1 = 0.1;
defaultOptions.w2 = 0.005;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = 0;
defaultOptions.TolX = 0;
defaultOptions.TolStagnationIteration = 30;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.mu_CR = [];
defaultOptions.initial.mu_F = [];
defaultOptions.initial.mu_H = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
R = options.R;
delta_CR = options.delta_CR;
delta_F = options.delta_F;
delta_H = options.delta_H;
w1 = options.w1;
w2 = options.w2;
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
	mu_CR = options.initial.mu_CR;
	mu_F = options.initial.mu_F;
	mu_H = options.initial.mu_H;
else
	X = [];
	f = [];
	mu_CR = [];
	mu_F = [];
	mu_H = [];
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
	'mu_F', 'mu_CR', 'mu_H');

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
	mu_F = options.F * ones(1, NP);
end

% mu_H
if isempty(mu_H)
	mu_H = options.H * ones(1, NP);
end

% mu_CR
if isempty(mu_CR)
	mu_CR = options.CR;
end

% Initialize variables
V = X(:, 1 : NP);
U = X(:, 1 : NP);
XT = X;
VT = X;
UT = X;
P = ones(NP, NP * NP);
r1r2 = zeros(1, NP);

% Display
if isDisplayIter
	displayitermessages(...
		X(:, 1 : NP), U(:, 1 : NP), f(1 : NP), countiter, XX, YY, ZZ, ...
		'mu_F', mean(mu_F), 'mu_CR', mu_CR, 'mu_H', mean(mu_H));
end

% Record
out = updateoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
	'mu_F', mean(mu_F), 'mu_CR', mu_CR, 'mu_H', mean(mu_H));

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
	S_CR = zeros(1, NP);
	CR = mu_CR + delta_CR * randn(1, NP);
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	
	F = zeros(1, NP);
	H = zeros(1, NP);
	
	for i = 1 : NP
		F(i) = mu_F(i) + mu_F(i) * delta_F * randn;
		if F(i) > 1
			F(i) = 1;
		elseif F(i) <= 0
			F(i) = 0.01 * mu_F(i) * (1 + rand);
		end
				
		H(i) = mu_H(i) + mu_H(i) * delta_H * randn;
		if H(i) > 1
			H(i) = 1;
		elseif H(i) <= 0
			H(i) = 0.01 * mu_H(i) * (1 + rand);
		end
	end
	
	A_Counter = 0;
	
	% Probability table
	
	for i = 1 : NP
		for j = 1 : NP
			P(i, j * NP + j) = 0;
		end
		
		P(i, :) = P(i, :) / sum(P(i, :));
	end
	
	% Mutation
	[~, besti] = min(f);
	for i = 1 : NP
		ptable = cumsum([0, P(i, :)]);
		for checkRetry = 1 : 3			
			% Generate r1, r2
			r1r2(i) = sum(rand >= ptable);
			r1 = 1 + floor((r1r2(i) - 1)/ NP);
			r2 = 1 + mod(r1r2(i) - 1, NP);
			
			V(:, i) = X(:, besti) + ...
				F(i) .* (X(:, i) - X(:, besti)) + ...
				H(i) .* (X(:, r1) - X(:, r2));
			
			% Check boundary
			if all(V(:, i) > lb) && all(V(:, i) < ub)
				break;
			end
		end
	end
	
	[B, ~] = eig(cov(X'));
	for i = 1 : NP
		if rand < R
			% Rotational Crossover
			XT(:, i) = B' * X(:, i);
			VT(:, i) = B' * V(:, i);
			jrand = floor(1 + D * rand);
			
			for j = 1 : D
				if rand < CR(i) || j == jrand
					UT(j, i) = VT(j, i);
				else
					UT(j, i) = XT(j, i);
				end
			end
			
			U(:, i) = B * UT(:, i);
		else
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
			'mu_F', mean(mu_F), 'mu_CR', mu_CR, 'mu_H', mean(mu_H));
	end
	
	% Selection
	FailedIteration = true;
	for i = 1 : NP
		fui = feval(fitfun, U(:, i));
		counteval = counteval + 1;
		
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
			S_CR(A_Counter + 1) = CR(i);
			A_Counter = A_Counter + 1;
			FailedIteration = false;
			P(i, r1r2(i)) = P(i, r1r2(i)) * 2;
			mu_F(i) = sqrt((1 - w2) * mu_F(i).^2 + w2 * F(i).^2);
			mu_H(i) = sqrt((1 - w2) * mu_H(i).^2 + w2 * H(i).^2);
		else
			P(i, r1r2(i)) = P(i, r1r2(i)) / 2;
			mu_F(i) = (1 - w2) * mu_F(i);
			mu_H(i) = (1 - w2) * mu_H(i);
		end
	end	
	
	% Update CR and F
	if A_Counter > 0
		mu_CR = (1 - w1) * mu_CR + w1 * mean(S_CR(1 : A_Counter));
	end
	
	% Record
	out = updateoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
		'mu_F', mean(mu_F), 'mu_CR', mu_CR, 'mu_H', mean(mu_H));
	
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

final.mu_F = mu_F;
final.mu_H = mu_H;
final.mu_CR = mu_CR;

out = finishoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
	'final', final, ...
	'mu_F', mean(mu_F), 'mu_CR', mu_CR, 'mu_H', mean(mu_H));
end
