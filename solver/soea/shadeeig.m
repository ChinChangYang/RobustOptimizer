function [xmin, fmin, out] = shadeeig(fitfun, lb, ub, maxfunevals, options)
% SHADEEIG SHADE/EIG algorithm
% SHADEEIG(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% SHADEEIG(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.H = 100;
defaultOptions.F = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.R = 0.5;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.MCR = [];
defaultOptions.initial.MF = [];

options = setdefoptions(options, defaultOptions);
NP = options.NP;
H = options.H;
R = options.R;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	fx = options.initial.f;
	A = options.initial.A;
	MCR = options.initial.MCR;
	MF = options.initial.MF;
else
	X = [];
	fx = [];
	A = [];
	MCR = [];
	MF = [];
end

D = numel(lb);
if ~isempty(X)
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

% Evaluation
if isempty(fx)
	fx = zeros(1, NP);
	for i = 1 : NP
		fx(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
end

% Sort
[fx, fidx] = sort(fx);
X = X(:, fidx);

% MF
if isempty(MF)
	MF = options.F * ones(H, 1);
end

% MCR
if isempty(MCR)
	MCR = options.CR * ones(H, 1);
end

% Initialize variables
V = X;
U = X;
XT = X;
VT = X;
UT = X;
C = cov(X');
k = 1;
r = zeros(1, NP);
p = zeros(1, NP);
pmin = 2 / NP;
A_size = 0;
fu = zeros(1, NP);
S_CR = zeros(1, NP);	% Set of crossover rate
S_F = zeros(1, NP);		% Set of scaling factor
S_df = zeros(1, NP);	% Set of df

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, ...
	'MF', mean(MF), ...
	'MCR', mean(MCR));

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(fx) <= ftarget;
	stagnation = countStagnation >= TolStagnationIteration;
	
	% Convergence conditions	
	if outofmaxfunevals || reachftarget || stagnation
		break;
	end
	
	% Reset S
	nS = 0;
	
	% Crossover rates
	CR = zeros(1, NP);	
	for i = 1 : NP
		r(i) = floor(1 + H * rand);
		CR(i) = MCR(r(i)) + 0.1 * randn;
	end
	
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	
	% Scaling factors
	F = zeros(1, NP);
	for i = 1 : NP
		while F(i) <= 0
			F(i) = cauchyrnd(MF(r(i)), 0.1);
		end
		
		if F(i) > 1
			F(i) = 1;
		end
	end
	
	% pbest
	for i = 1 : NP
		p(i) = pmin + rand * (0.2 - pmin);
	end
	
	XA = [X, A];
	
	% Mutation
	for i = 1 : NP
		% Generate pbest_idx
		pbest = floor(1 + round(p(i) * NP) * rand);
		
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
		
		V(:, i) = X(:, i) + F(i) .* (X(:, pbest) - X(:, i)) ...
			+ F(i) .* (X(:, r1) - XA(:, r2));
	end
	
	[B, ~] = eig(C);
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
			nS = nS + 1;
			S_CR(nS)	= CR(i);
			S_F(nS)		= F(i);
			S_df(nS)	= abs(fu(i) - fx(i));
			X(:, i)		= U(:, i);
			fx(i)		= fu(i);
			
			if A_size < NP
				A = [A, X(:, i)];
				A_size = A_size + 1;
			else
				ri = floor(1 + NP * rand);
				A(:, ri) = X(:, i);
			end
						
			FailedIteration = false;
		elseif fu(i) == fx(i)
			X(:, i) = U(:, i);
			fx(i)	= fu(i);
		end
	end
	
	% Update MCR and MF
	if nS > 0
		w = S_df(1 : nS) ./ sum(S_df(1 : nS));
		MCR(k) = sum(w .* S_CR(1 : nS));
		MF(k) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));
		k = k + 1;
		if k > H
			k = 1;
		end
	end
	
	% Update C
	C = cov(X');
	
	% Sort	
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	
	% Record
	out = updateoutput(out, X, fx, counteval, ...
		'MF', mean(MF), ...
		'MCR', mean(MCR));
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end	
end

fmin = fx(1);
xmin = X(:, 1);

if fmin < out.bestever.fmin
	out.bestever.fmin = fmin;
	out.bestever.xmin = xmin;
end

final.A = A;
final.MCR = MCR;
final.MF = MF;

out = finishoutput(out, X, fx, counteval, ...
	'final', final, ...
	'MF', mean(MF), ...
	'MCR', mean(MCR));
end
