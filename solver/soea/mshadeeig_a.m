function [xmin, fmin, out] = mshadeeig_a(fitfun, lb, ub, maxfunevals, options)
% MSHADEEIG_A Mutable SHADE/EIG algorithm (Alternative) 
% MSHADEEIG_A(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% MSHADEEIG_A(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 114;
defaultOptions.H = 114;
defaultOptions.F = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.R = 0.5;
defaultOptions.cc = 0.1;
defaultOptions.pmin = 2/100;
defaultOptions.pmax = 0.2;
defaultOptions.Q = 50;
defaultOptions.deltaF = 0.1;
defaultOptions.deltaCR = 0.1;
defaultOptions.deltaR = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = 0;
defaultOptions.TolX = 0;
defaultOptions.TolStagnationIteration = 100;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.MCR = [];
defaultOptions.initial.MF = [];
defaultOptions.initial.MR = [];

options = setdefoptions(options, defaultOptions);
NP = options.NP;
H = options.H;
pmin = options.pmin;
pmax = options.pmax;
cc = options.cc;
Q = options.Q;
deltaF = options.deltaF;
deltaCR = options.deltaCR;
deltaR = options.deltaR;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	fx = options.initial.f;
	A = options.initial.A;
	MCR = options.initial.MCR;
	MF = options.initial.MF;
	MR = options.initial.MR;
else
	X = [];
	fx = [];
	A = [];
	MCR = [];
	MF = [];
	MR = [];
end

D = numel(lb);
if ~isempty(X)
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'MF', 'MCR', 'MR', 'FC1Q', 'FCMEDIAN', 'FC3Q');

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

% MF
if isempty(MF)
	MF = options.F * ones(H, 1);
end

% MCR
if isempty(MCR)
	MCR = options.CR * ones(H, 1);
end

% MR
if isempty(MR)
	MR = options.R * ones(H, 1);
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
rf = zeros(1, NP);
A_size = 0;
fu = zeros(1, NP);
S_CR = zeros(1, NP);	% Set of crossover rates
S_F = zeros(1, NP);		% Set of scaling factors
S_df = zeros(1, NP);	% Set of df
S_R = zeros(1, NP);		% Set of eigenvector ratio
FC = zeros(1, NP);		% Fail Counter

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, ...
	'MF', mean(MF), ...
	'MCR', mean(MCR), ...
	'MR', mean(MR), ...
	'FC1Q', quantile(FC, 0.25), ...
	'FCMEDIAN', median(FC), ...
	'FC3Q', quantile(FC, 0.75));

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(fx) <= ftarget;
	
	% Convergence conditions	
	if outofmaxfunevals
		out.stopflag = 'outofmaxfunevals';
		break;
	elseif reachftarget
		out.stopflag = 'reachftarget';
		break;
	end
	
	% Reset S
	nS = 0;
	
	% Crossover rates
	CR = zeros(1, NP);	
	for i = 1 : NP
		r(i) = floor(1 + H * rand);
		CR(i) = MCR(r(i)) + deltaCR * randn;
	end
	
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	
	% Scaling factors
	F = zeros(1, NP);
	for i = 1 : NP
		while F(i) <= 0
			F(i) = cauchyrnd(MF(r(i)), deltaF);
		end
		
		if F(i) > 1
			F(i) = 1;
		end
	end
	
	% Eigenvector ratio
	R = zeros(1, NP);	
	for i = 1 : NP
		R(i) = MR(r(i)) + deltaR * randn;
	end
	
	R(R > 1) = 1;
	R(R < 0) = 0;
	
	% pbest
	for i = 1 : NP
		p(i) = pmin + rand * (pmax - pmin);
	end
	
	% Archive
	XA = [X, A];
	
	% Emergency mode
	for i = 1 : NP		
		if FC(i) <= Q
			rf(i) = i;
		else
			rf(i) = floor(1 + NP * rand);
		end
	end
	
	% Mutation
	for i = 1 : NP
		% Generate pbest_idx
		pbest = floor(1 + round(p(i) * NP) * rand);
		
		% Generate r1
		r1 = floor(1 + NP * rand);
		while rf(i) == r1
			r1 = floor(1 + NP * rand);
		end
		
		% Generate r2
		r2 = floor(1 + (NP + A_size) * rand);
		while rf(i) == r1 || r1 == r2
			r2 = floor(1 + (NP + A_size) * rand);
		end
				
		V(:, i) = X(:, rf(i)) + F(rf(i)) .* (X(:, pbest) - X(:, rf(i))) ...
			+ F(rf(i)) .* (X(:, r1) - XA(:, r2));
	end
	
	[B, ~] = eig(C);
	for i = 1 : NP
		if rand < R(rf(i))
			% Rotational Crossover
			XT(:, i) = B' * X(:, rf(i));
			VT(:, i) = B' * V(:, i);
			jrand = floor(1 + D * rand);			
			for j = 1 : D
				if rand < CR(rf(i)) || j == jrand
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
				if rand < CR(rf(i)) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, rf(i));
				end
			end
		end
	end
	
	% Correction for outside of boundaries
	for i = 1 : NP
		for j = 1 : D
			if U(j, i) < lb(j)
				U(j, i) = 0.5 * (lb(j) + X(j, rf(i)));
			elseif U(j, i) > ub(j)
				U(j, i) = 0.5 * (ub(j) + X(j, rf(i)));
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
	for i = 1 : NP		
		if fu(i) < fx(i)
			nS = nS + 1;
			S_CR(nS)	= CR(rf(i));
			S_F(nS)		= F(rf(i));
			S_df(nS)	= abs(fu(i) - fx(i));
			S_R(nS)		= R(rf(i));
			X(:, i)		= U(:, i);
			fx(i)		= fu(i);
			FC(i)		= 0;
			
			if A_size < NP
				A_size = A_size + 1;
				A(:, A_size) = X(:, i);
			else
				ri = floor(1 + NP * rand);
				A(:, ri) = X(:, i);
			end
		else
			FC(i) = FC(i) + 1;
		end
	end
	
	% Update MCR and MF
	if nS > 0
		w = S_df(1 : nS) ./ sum(S_df(1 : nS));
		MCR(k) = sum(w .* S_CR(1 : nS));
		MF(k) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));		
		MR(k) = sum(w .* S_R(1 : nS));
		k = k + 1;
		if k > H
			k = 1;
		end
	end
	
	% Update C
	C = (1 - cc) * C + cc * cov(X');
	
	% Sort	
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	FC = FC(fidx);
	
	% Record
	out = updateoutput(out, X, fx, counteval, ...
		'MF', mean(MF), ...
		'MCR', mean(MCR), ...
		'MR', mean(MR), ...
		'FC1Q', quantile(FC, 0.25), ...
		'FCMEDIAN', median(FC), ...
		'FC3Q', quantile(FC, 0.75));
	
	% Iteration counter
	countiter = countiter + 1;
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
final.MR = MR;

out = finishoutput(out, X, fx, counteval, ...
	'final', final, ...
	'MF', mean(MF), ...
	'MCR', mean(MCR), ...
	'MR', mean(MR), ...
	'FC1Q', quantile(FC, 0.25), ...
	'FCMEDIAN', median(FC), ...
	'FC3Q', quantile(FC, 0.75));
end
