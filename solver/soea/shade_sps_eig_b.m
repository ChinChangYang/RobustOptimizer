function [xmin, fmin, out] = shade_sps_eig_b(fitfun, lb, ub, maxfunevals, options)
% SHADE_SPS_EIG SHADE algorithm with SPS and EIG Framework (B)
% SHADE_SPS_EIG(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% SHADE_SPS_EIG(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.Q = 70;
defaultOptions.R = 0.5;
defaultOptions.wc = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.MCR = [];
defaultOptions.initial.MF = [];
defaultOptions.initial.MR = [];
defaultOptions.ConstraintHandling = 'Interpolation';

options = setdefoptions(options, defaultOptions);
Q = options.Q;
wc = options.wc;
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
if isempty(X)
	NP = options.NP;
	H = NP;
else
	[~, NP] = size(X);
	H = NP;
end

% Initialize variables
counteval = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'FC', ...
	'muMF', ...
	'muMCR', ...
	'muMR');

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
if isempty(MF)
	MF = options.F * ones(H, 1);
end

% mu_CR
if isempty(MCR)
	MCR = options.CR * ones(H, 1);
end

% mu_R
if isempty(MR)
	MR = options.R * ones(H, 1);
end

% Initialize variables
V = X;
U = X;
XT = X;
VT = X;
UT = X;
SPT = X;
C = cov(X');
k = 1;
r = zeros(1, NP);
p = zeros(1, NP);
pmin = 2 / NP;
A_size = 0;
fu = zeros(1, NP);
S_CR = zeros(1, NP);	% Set of crossover rate
S_F = zeros(1, NP);		% Set of scaling factor
S_R = zeros(1, NP);		% Set of EIG ratio
S_df = zeros(1, NP);	% Set of df
FC = zeros(1, NP);		% Consecutive Failure Counter
Chy = cauchyrnd(0, 0.1, NP + 10);
iChy = 1;
SP = X;
fSP = fx;
iSP = 1;
[~, sortidxfSP] = sort(fSP);

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter, ...
	'FC', FC, ...
	'muMF', mean(MF), ...
	'muMCR', mean(MCR), ...
	'muMR', mean(MR));

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
			F(i) = MF(r(i)) + Chy(iChy);
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
	
	% EIG ratio
	R = zeros(1, NP);
	
	for i = 1 : NP
		r(i) = floor(1 + H * rand);
		R(i) = MR(r(i)) + 0.1 * randn;
	end
	
	R(R > 1) = 1;
	R(R < 0) = 0;	
	
	% pbest
	for i = 1 : NP
		p(i) = pmin + rand * (0.2 - pmin);
	end
	
	XA = [X, A];
	SPA = [SP, A];
	
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
		
		if FC(i) <= Q
			V(:, i) = X(:, i) + F(i) .* (X(:, pbest) - X(:, i)) ...
				+ F(i) .* (X(:, r1) - XA(:, r2));
		else
			V(:, i) = SP(:, i) + ...
				F(i) .* (SP(:, sortidxfSP(pbest)) - SP(:, i)) ...
				+ F(i) .* (SP(:, r1) - SPA(:, r2));
		end
	end
	
	C = (1 - wc) * C + wc * cov(X');
	[B, ~] = eig(C);
	for i = 1 : NP
		jrand = floor(1 + D * rand);
		
		if rand < R(i)
			if FC(i) <= Q
				% EIG Framework
				XT(:, i) = B' * X(:, i);
				VT(:, i) = B' * V(:, i);
				for j = 1 : D
					if rand < CR(i) || j == jrand
						UT(j, i) = VT(j, i);
					else
						UT(j, i) = XT(j, i);
					end
				end
			else
				% EIG+SPS Framework
				XT(:, i) = B' * X(:, i);
				SPT(:, i) = B' * SP(:, i);
				for j = 1 : D
					if rand < CR(i) || j == jrand
						UT(j, i) = VT(j, i);
					else
						UT(j, i) = SPT(j, i);
					end
				end
			end				
			U(:, i) = B * UT(:, i);
		else			
			if FC(i) <= Q
				% Binomial Crossover
				for j = 1 : D
					if rand < CR(i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = X(j, i);
					end
				end
			else
				% SPS Framework
				for j = 1 : D
					if rand < CR(i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = SP(j, i);
					end
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
			nS = nS + 1;
			S_CR(nS)	= CR(i);
			S_F(nS)		= F(i);
			S_R(nS)		= R(i);
			S_df(nS)	= abs(fu(i) - fx(i));
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
	
	% Sort
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	FC = FC(fidx);
	[~, sortidxfSP] = sort(fSP);
	
	% Record
	out = updateoutput(out, X, fx, counteval, countiter, ...
		'FC', FC, ...
		'muMF', mean(MF), ...
		'muMCR', mean(MCR), ...
		'muMR', mean(MR));
	
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
	'FC', zeros(NP, 1), ...
	'muMF', mean(MF), ...
	'muMCR', mean(MCR), ...
	'muMR', mean(MR));
end
