function [xmin, fmin, out] = jade_sps(fitfun, lb, ub, maxfunevals, options)
% JADE_SPS JADE algorithm with SPS Framework
% JADE_SPS(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% JADE_SPS(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.Q = 70;
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
defaultOptions.ConstraintHandling = 'Interpolation';

options = setdefoptions(options, defaultOptions);
Q = options.Q;
delta_CR = options.delta_CR;
delta_F = options.delta_F;
p = options.p;
w = options.w;
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
	A = options.initial.A;
	mu_CR = options.initial.mu_CR;
	mu_F = options.initial.mu_F;
else
	X = [];
	fx = [];
	A = [];
	mu_CR = [];
	mu_F = [];
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
	'mu_F', 'mu_CR', ...
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

% Initialize archive
if isempty(A)
	A = X;
end

% Sort
[fx, fidx] = sort(fx);
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
fu = zeros(1, NP);
S_CR = zeros(1, NP);	% Set of crossover rate
S_F = zeros(1, NP);		% Set of scaling factor
FC = zeros(1, NP);		% Consecutive Failure Counter
SP = X;
fSP = fx;
iSP = 1;
[~, sortidxfSP] = sort(fSP);
Chy = cauchyrnd(0, delta_F, NP + 10);
iChy = 1;

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter, ...
	'mu_F', mu_F, ...
	'mu_CR', mu_CR, ...
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
	
	% Reset S
	nS = 0;
		
	% Crossover rates
	CR = mu_CR + delta_CR * randn(1, NP);
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	
	% Scaling factors
	F = zeros(1, NP);
	for i = 1 : NP
		while F(i) <= 0
			F(i) = mu_F + Chy(iChy);
			if iChy < numel(Chy)
				iChy = iChy + 1;
			else
				iChy = 1;
			end
		end
		
		if F(i) > 1
			F(i) = 1;
		end
	end
	
	XA = [X, A];
	SPA = [SP, A];
	
	% Mutation
	for i = 1 : NP				
		% Generate pbest_idx
		pbest = max(1, ceil(rand * pbest_size));
		
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
			
		if FC(i) <= Q
			V(:, i) = X(:, i) + F(i) .* (X(:, pbest) - X(:, i)) ...
				+ F(i) .* (X(:, r1) - XA(:, r2));
		else
			V(:, i) = SP(:, i) + F(i) .* (SP(:, sortidxfSP(pbest)) - SP(:, i)) ...
				+ F(i) .* (SP(:, r1) - SPA(:, r2));
		end
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		if FC(i) <= Q
			for j = 1 : D
				if rand < CR(i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, i);
				end
			end
		else
			for j = 1 : D
				if rand < CR(i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = SP(j, i);
				end
			end			
		end
	end
	
	if interpolation
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
		if fu(i) < fx(i)
			nS			= nS + 1;
			S_CR(nS)	= CR(i);
			S_F(nS)		= F(i);
			X(:, i)		= U(:, i);
			fx(i)		= fu(i);
			SP(:, iSP)	= U(:, i);
			fSP(iSP)	= fu(i);
			iSP			= mod(iSP, NP) + 1;
			FC(i)		= 0;
			
			if A_size < NP
				A_size = A_size + 1;
				A(:, A_size) = X(:, i);
			else
				ri = floor(1 + NP * rand);
				A(:, ri) = X(:, i);
			end
			
			FailedIteration = false;
		else
			FC(i) = FC(i) + 1;
		end
	end
	
	% Update CR and F
	if nS > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : nS));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : nS).^2) / sum(S_F(1 : nS));
	end
	
	% Sort	
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	FC = FC(fidx);
	[~, sortidxfSP] = sort(fSP);
	
	% Record
	out = updateoutput(out, X, fx, counteval, countiter, ...
		'mu_F', mu_F, ...
		'mu_CR', mu_CR, ...
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

final.A = A;
final.mu_F = mu_F;
final.mu_CR = mu_CR;

out = finishoutput(out, X, fx, counteval, countiter, ...
	'final', final, ...
	'mu_F', mu_F, ...
	'mu_CR', mu_CR, ...
	'FC', zeros(NP, 1));
end
