function [xmin, fmin, out] = lshade_sps(fitfun, lb, ub, maxfunevals, options)
% LSHADE_SPS L-SHADE algorithm with SPS framework
% LSHADE_SPS(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% LSHADE_SPS(..., options) minimize the function by solver options.
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
defaultOptions.Q = 64;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.usefunevals = inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.nA = [];
defaultOptions.initial.MCR = [];
defaultOptions.initial.MF = [];
defaultOptions.initial.iM = [];
defaultOptions.initial.FC = [];
defaultOptions.initial.SP = [];
defaultOptions.initial.fSP = [];
defaultOptions.initial.iSP = [];
defaultOptions.initial.counteval = [];
defaultOptions.initial.countiter = [];
defaultOptions.initial.countstagnation = [];
defaultOptions.ConstraintHandling = 'Interpolation';
defaultOptions.EarlyStop = 'fitness';

options = setdefoptions(options, defaultOptions);
Ar = options.Ar;
p = options.p;
H = options.H;
NPmin = eval(options.NPmin);
Q = options.Q;
NPinit = options.NP;
usefunevals = options.usefunevals;
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
	X		= options.initial.X;
	fx		= options.initial.f;
	A		= options.initial.A;
	nA		= options.initial.nA;
	MCR		= options.initial.MCR;
	MF		= options.initial.MF;
	iM		= options.initial.iM;
	FC		= options.initial.FC;
	SP		= options.initial.SP;
	fSP		= options.initial.fSP;
	iSP		= options.initial.iSP;
	counteval = options.initial.counteval;
	countiter = options.initial.countiter;
	countstagnation = options.initial.countstagnation;
else
	X		= [];
	fx		= [];
	A		= [];
	nA		= [];
	MCR		= [];
	MF		= [];
	iM		= [];
	FC		= [];
	SP		= [];
	fSP		= [];
	iSP		= [];
	counteval = [];
	countiter = [];
	countstagnation = [];
end

D = numel(lb);
if isempty(X)
	NP = options.NP;
else
	[~, NP] = size(X);
end

% Initialize out
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'muMF', ...
	'muMCR', ...
	'muFC');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
end

% counteval
if isempty(counteval)
	counteval = 0;
end

% countiter
if isempty(countiter)
	countiter = 1;
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
	Asize = round(Ar * NP);
	A = zeros(D, Asize);
	for i = 1 : Asize
		A(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
else
	[~, Asize] = size(A);
	if Asize > round(Ar * NP)
		Asize = round(Ar * NP);
		A = A(:, 1 : Asize);
	elseif Asize < round(Ar * NP)
		Asize = round(Ar * NP);
		A = zeros(D, Asize);
		for i = 1 : Asize
			A(:, i) = lb + (ub - lb) .* rand(D, 1);
		end
	end
end

if isempty(nA)
	nA = 0;
else
	nA = min(nA, Asize);
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

% iM
if isempty(iM)
	iM = 1;
end

% FC
if isempty(FC)
	FC = zeros(1, NP);		% Consecutive Failure Counter
end

% SP and fSP
if isempty(SP)
	SP = X;
	fSP = fx;
elseif isempty(fSP)
	fSP = zeros(1, NP);
	for i = 1 : NP
		fSP(i) = feval(fitfun, SP(:, i));
		counteval = counteval + 1;
	end
end

% iSP
if isempty(iSP)
	iSP = 1;
end

% Initialize variables
V = X;
U = X;
fu = zeros(1, NP);
S_CR = zeros(1, NP);	% Set of crossover rate
S_F = zeros(1, NP);		% Set of scaling factor
S_df = zeros(1, NP);	% Set of df
Chy = cauchyrnd(0, 0.1, NP + 10);
iChy = 1;
[~, sortidxfSP] = sort(fSP);

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter, ...
	'muMF', mean(MF), ...
	'muMCR', mean(MCR), ...
	'muFC', mean(FC));

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	outofusefunevals = counteval > usefunevals - NP;
	if ~EarlyStopOnFitness && ~AutoEarlyStop
		if outofmaxfunevals || outofusefunevals
			break;
		end
	elseif AutoEarlyStop
		reachftarget = min(fx) <= ftarget;
		TolX = 10 * eps(mean(X(:)));
		solutionconvergence = std(X(:)) <= TolX;
		TolFun = 10 * eps(mean(fx));
		functionvalueconvergence = std(fx(:)) <= TolFun;
		stagnation = countstagnation >= TolStagnationIteration;
		
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
	
	% Iteration counter
	countiter = countiter + 1;
	
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
	end
	F(F > 1) = 1;
	
	% pbest
	pbest = 1 + floor(max(2, round(p * NP)) * rand(1, NP));
	
	% Population + archive
	XA = [X, A];
	SPA = [SP, A];
	
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
		
		if FC(i) <= Q
			V(:, i) = X(:, i) ...
				+ F(i) .* (X(:, pbest(i)) - X(:, i)) ...
				+ F(i) .* (X(:, r1) - XA(:, r2));
		else
			V(:, i) = SP(:, i) ...
				+ F(i) .* (SP(:, sortidxfSP(pbest(i))) - SP(:, i)) ...
				+ F(i) .* (SP(:, r1) - SPA(:, r2));
		end
	end
	
	if interpolation
		% Correction for outside of boundaries
		for i = 1 : NP
			for j = 1 : D
				if V(j, i) < lb(j)
					if FC(i) <= Q
						V(j, i) = 0.5 * (lb(j) + X(j, i));
					else
						V(j, i) = 0.5 * (lb(j) + SP(j, i));
					end
				elseif V(j, i) > ub(j)
					if FC(i) <= Q
						V(j, i) = 0.5 * (ub(j) + X(j, i));
					else
						V(j, i) = 0.5 * (ub(j) + SP(j, i));
					end
				end
			end
		end
	end
	
	for i = 1 : NP
		jrand = floor(1 + D * rand);
		
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
			FC(i)		= 0;		
			SP(:, iSP)	= U(:, i);
			fSP(iSP)	= fu(i);
			iSP			= mod(iSP, NP) + 1;	
		elseif fu(i) == fx(i)			
			X(:, i)		= U(:, i);
			FC(i)		= FC(i) + 1;
		else			
			FC(i)		= FC(i) + 1;
		end
	end
	
	% Update MCR and MF
	if nS > 0
		w = S_df(1 : nS) ./ sum(S_df(1 : nS));
		
		if all(S_CR(1 : nS) == 0)
			MCR(iM) = -1;
		elseif MCR(iM) ~= -1
			MCR(iM) = sum(w .* S_CR(1 : nS) .* S_CR(1 : nS)) / sum(w .* S_CR(1 : nS));
		end
		
		MF(iM) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));
		iM = mod(iM, H) + 1;
	end
	
	% Sort
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	FC = FC(fidx);
	
	% Update NP and population
    NP = round(NPinit - (NPinit - NPmin) * counteval / maxfunevals);
    fx = fx(1 : NP);
    X = X(:, 1 : NP);
    Asize = round(Ar * NP);	
	if nA > Asize
		nA = Asize;
		A = A(:, 1 : Asize);
	end
    FC = FC(1 : NP);
	[~, sortidxfSP] = sort(fSP);    
    remainingfSPidx = sortidxfSP <= NP;
    SP = SP(:, remainingfSPidx);
    fSP = fSP(:, remainingfSPidx);
    sortidxfSP = sortidxfSP(remainingfSPidx);
    iSP	= mod(iSP - 1, NP) + 1;
	
	% Record
	out = updateoutput(out, X, fx, counteval, countiter, ...
		'muMF', mean(MF), ...
		'muMCR', mean(MCR), ...
		'muFC', mean(FC));
	
	% Stagnation iteration
	if FailedIteration
		countstagnation = countstagnation + 1;
	else
		countstagnation = 0;
	end
end

fmin = fx(1);
xmin = X(:, 1);

final.A			= A;
final.nA		= nA;
final.MCR		= MCR;
final.MF		= MF;
final.iM		= iM;
final.FC		= FC;
final.SP		= SP;
final.fSP		= fSP;
final.iSP		= iSP;
final.counteval = counteval;
final.countiter = countiter;
final.countstagnation = countstagnation;

out = finishoutput(out, X, fx, counteval, countiter, ...
	'final', final, ...
	'muMF', mean(MF), ...
	'muMCR', mean(MCR), ...
	'muFC', mean(FC));
end
