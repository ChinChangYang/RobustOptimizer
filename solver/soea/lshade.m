function [xmin, fmin, out] = lshade(fitfun, lb, ub, maxfunevals, options)
% LSHADE L-SHADE algorithm
% LSHADE(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% LSHADE(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 540;
defaultOptions.F = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.Ar = 2.6;
defaultOptions.p = 0.11;
defaultOptions.H = 6;
defaultOptions.NPmin = '4';
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.MCR = [];
defaultOptions.initial.MF = [];
defaultOptions.ConstraintHandling = 'Interpolation';
defaultOptions.EarlyStop = 'fitness';

options = setdefoptions(options, defaultOptions);
Ar = options.Ar;
p = options.p;
H = options.H;
NPmin = eval(options.NPmin);
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

if isequal(options.ConstraintHandling, 'Interpolation')
	interpolation = true;
else
	interpolation = false;
end


if ~isempty(strfind(options.EarlyStop, 'fitness'))
	EarlyStopOnFitness = true;
	AutoEarlyStop = false;
elseif ~isempty(strfind(options.EarlyStop, 'auto'))
	EarlyStopOnFitness = false;
	AutoEarlyStop = true;
else
	EarlyStopOnFitness = false;
	AutoEarlyStop = false;
end

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
	'muMF', ...
	'muMCR');

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

% Archive
if isempty(A)
	Asize = round(Ar * NP);
	A = zeros(D, Asize);
	for i = 1 : Asize
		A(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
	nA = 0;
else
	[~, Asize] = size(A);
	Ar = Asize / NP;
	nA = Asize;
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

% Initialize variables
V = X;
U = X;
k = 1;
fu = zeros(1, NP);
S_CR = zeros(1, NP);	% Set of crossover rate
S_F = zeros(1, NP);		% Set of scaling factor
S_df = zeros(1, NP);	% Set of df
Chy = cauchyrnd(0, 0.1, NP + 10);
iChy = 1;
NPinit = NP;

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter, ...
	'muMF', mean(MF), ...
	'muMCR', mean(MCR));

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	if ~EarlyStopOnFitness && ~AutoEarlyStop
		if outofmaxfunevals
			break;
		end
	elseif AutoEarlyStop
		reachftarget = min(fx) <= ftarget;
		TolX = 10 * eps(mean(X(:)));
		solutionconvergence = std(X(:)) <= TolX;
		TolFun = 10 * eps(mean(fx));
		functionvalueconvergence = std(fx(:)) <= TolFun;
		stagnation = countStagnation >= TolStagnationIteration;
		
		if outofmaxfunevals || ...
				reachftarget || ...
				solutionconvergence || ...
				functionvalueconvergence || ...
				stagnation
			break;
		end
	elseif EarlyStopOnFitness
		reachftarget = min(fx) <= ftarget;
		
		if outofmaxfunevals || ...
				reachftarget
			break;
		end
	end
	
	% Memory Indices
	r = floor(1 + H * rand(1, NP));
	
	% Crossover rates
	CR = MCR(r)' + 0.1 * randn(1, NP);	
	CR((CR < 0) | (MCR(r)' == -1)) = 0;
	CR(CR > 1) = 1;
	
	% Scaling factors
	F = zeros(1, NP);
	for i = 1 : NP
		while F(i) <= 0
			F(i) = MF(r(i)) + Chy(iChy);
			iChy = mod(iChy, numel(Chy)) + 1;
		end
		
		if F(i) > 1
			F(i) = 1;
		end
	end
	
	% pbest
	pbest = 1 + floor(max(2, round(p * NP)) * rand(1, NP));
	
	% Population + archive
	XA = [X, A];
	
	% Mutation
	for i = 1 : NP		
		% Generate r1
		r1 = floor(1 + NP * rand);
		while i == r1
			r1 = floor(1 + NP * rand);
		end
		
		% Generate r2
		r2 = floor(1 + (NP + nA) * rand);
		while i == r1 || r1 == r2
			r2 = floor(1 + (NP + nA) * rand);
		end
		
		V(:, i) = X(:, i) + F(i) .* (X(:, pbest(i)) - X(:, i)) ...
			+ F(i) .* (X(:, r1) - XA(:, r2));
	end
	
	if interpolation
		% Correction for outside of boundaries
		for i = 1 : NP
			for j = 1 : D
				if V(j, i) < lb(j)
					V(j, i) = 0.5 * (lb(j) + X(j, i));
				elseif V(j, i) > ub(j)
					V(j, i) = 0.5 * (ub(j) + X(j, i));
				end
			end
		end
	end
	
	for i = 1 : NP
		jrand = floor(1 + D * rand);
		
		% Binomial Crossover
		for j = 1 : D
			if rand < CR(i) || j == jrand
				U(j, i) = V(j, i);
			else
				U(j, i) = X(j, i);
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
	nS = 0;
	for i = 1 : NP
		if fu(i) < fx(i)
			nS			= nS + 1;
			S_CR(nS)	= CR(i);
			S_F(nS)		= F(i);
			S_df(nS)	= abs(fu(i) - fx(i));
			X(:, i)		= U(:, i);
			fx(i)		= fu(i);
			
			if nA < Asize
				A(:, nA + 1)	= X(:, i);
				nA				= nA + 1;
			else
				ri				= floor(1 + Asize * rand);
				A(:, ri)		= X(:, i);
			end
			
			FailedIteration = false;
		elseif fu(i) == fx(i)			
			X(:, i)		= U(:, i);
		end
	end
	
	% Update MCR and MF
	if nS > 0
		w = S_df(1 : nS) ./ sum(S_df(1 : nS));
		
		if all(S_CR(1 : nS) == 0)
			MCR(k) = -1;
		elseif MCR(k) ~= -1
			MCR(k) = sum(w .* S_CR(1 : nS) .* S_CR(1 : nS)) / sum(w .* S_CR(1 : nS));
		end
		
		MF(k) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));
		k = mod(k, H) + 1;
	end
	
	% Sort
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	
	% Update NP and population
    NP = round(NPinit - (NPinit - NPmin) * counteval / maxfunevals);
    fx = fx(1 : NP);
    X = X(:, 1 : NP);
    Asize = round(Ar * NP);	
	if nA > Asize
		nA = Asize;
	end
	
	% Record
	out = updateoutput(out, X, fx, counteval, countiter, ...
		'muMF', mean(MF), ...
		'muMCR', mean(MCR));
	
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
	'muMF', mean(MF), ...
	'muMCR', mean(MCR));
end
