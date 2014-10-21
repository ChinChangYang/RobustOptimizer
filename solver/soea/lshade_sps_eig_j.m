function [xmin, fmin, out] = lshade_sps_eig_j(fitfun, lb, ub, maxfunevals, options)
% LSHADE_SPS_EIG_J L-SHADE algorithm with SPS+EIG framework (variant J)
% LSHADE_SPS_EIG_J(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% LSHADE_SPS_EIG_J(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 540;
defaultOptions.F = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.ER = 1.0;
defaultOptions.p = 0.11;
defaultOptions.H = 6;
defaultOptions.Q = 64;
defaultOptions.Ar = 2.6;
defaultOptions.cw = 0.3;
defaultOptions.erw = 0.2;
defaultOptions.CRmin = 0.0;
defaultOptions.CRmax = 0.5;
defaultOptions.NPmin = '4';
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.usefunevals = inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.nA = [];
defaultOptions.initial.MER = [];
defaultOptions.initial.MCR = [];
defaultOptions.initial.MF = [];
defaultOptions.initial.psai = [];
defaultOptions.initial.iM = [];
defaultOptions.initial.FC = [];
defaultOptions.initial.C = [];
defaultOptions.initial.SP = [];
defaultOptions.initial.fSP = [];
defaultOptions.initial.iSP = [];
defaultOptions.initial.counteval = [];
defaultOptions.initial.countiter = [];
defaultOptions.initial.countstagnation = [];
defaultOptions.initial.countcon = [];
defaultOptions.ConstraintHandling = 'Interpolation';
defaultOptions.EpsilonValue = 0;
defaultOptions.nonlcon = [];
defaultOptions.EarlyStop = 'none';

options = setdefoptions(options, defaultOptions);
p = options.p;
H = options.H;
Q = options.Q;
cw = options.cw;
Ar = options.Ar;
NPinit = options.NP;
cwinit = options.cw;
erw = options.erw;
CRmin = options.CRmin;
CRmax = options.CRmax;
NPmin = eval(options.NPmin);
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

nonlcon = options.nonlcon;
EpsilonValue = options.EpsilonValue;
if ~isempty(strfind(options.ConstraintHandling, 'EpsilonMethod'))
	EpsilonMethod = true;
else
	EpsilonMethod = false;
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
	MER		= options.initial.MER;
	MF		= options.initial.MF;
	MCR		= options.initial.MCR;
	psai_x	= options.initial.psai;
	iM		= options.initial.iM;
	FC		= options.initial.FC;
	C		= options.initial.C;
	SP		= options.initial.SP;
	fSP		= options.initial.fSP;
	iSP		= options.initial.iSP;
	counteval = options.initial.counteval;
	countiter = options.initial.countiter;
	countstagnation = options.initial.countstagnation;
	countcon = options.initial.countcon;
else
	X		= [];
	fx		= [];
	A		= [];
	nA		= [];
	MER		= [];
	MF		= [];
	MCR		= [];
	psai_x	= [];
	iM		= [];
	FC		= [];
	C		= [];
	SP		= [];
	fSP		= [];
	iSP		= [];
	counteval = [];
	countiter = [];
	countstagnation = [];
	countcon = [];
end

D = numel(lb);
if isempty(X)
	NP = options.NP;
else
	[~, NP] = size(X);
end

% Initialize variables
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'countcon', ...
	'muMF', ...
	'muMCR', ...
	'muMER', ...
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

% countstagnation
if isempty(countstagnation)
	countstagnation = 0;
end

% countcon
if isempty(countcon)
	countcon = 0;
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

% Constraint violation
if isempty(psai_x) && EpsilonMethod
	psai_x = zeros(1, NP);
	for i = 1 : NP		
		clbx = lb - X(:, i);
		cubx = X(:, i) - ub;
		psai_x(i) = sum(clbx(clbx > 0)) + sum(cubx(cubx > 0));
		
		if ~isempty(nonlcon)			
			[cx, ceqx] = feval(nonlcon, X(:, i));
			countcon = countcon + 1;
			psai_x(i) = psai_x(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
		end
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
if ~EpsilonMethod
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
else
	PsaiFx = [psai_x', fx'];
	[~, SortingIndex] = sortrows(PsaiFx);
	X = X(:, SortingIndex);
	fx = fx(SortingIndex);
	psai_x = psai_x(SortingIndex);
end

% MF
if isempty(MF)
	MF = options.F * ones(H, 1);
end

% MCR
if isempty(MCR)
	MCR = options.CR * ones(H, 1);
end

% MER
if isempty(MER)
	MER = options.ER * ones(H, 1);
end

% iM
if isempty(iM)
	iM = 1;
end

% FC
if isempty(FC)
	FC = zeros(1, NP);		% Consecutive Failure Counter
end

% C
if isempty(C)
	C = cov(X');
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
S_ER = zeros(1, NP);	% Set of EIG rate
S_F = zeros(1, NP);		% Set of scaling factor
S_CR = zeros(1, NP);	% Set of crossover rate
S_df = zeros(1, NP);	% Set of df
Chy = cauchyrnd(0, 0.1, NP + 10);
iChy = 1;
psai_u = zeros(1, NP);
[~, sortidxfSP] = sort(fSP);

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter, ...
	'countcon', countcon, ...
	'muMF', mean(MF), ...
	'muMCR', mean(MCR), ...
	'muMER', mean(MER), ...
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
	
	% EIG rates
	ER = MER(r)' + erw * randn(1, NP);
	
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
	
	% Crossover rates	
	CR = rand(1, NP);
	CR(CR <= MCR(r)') = CRmin;
	CR(CR > MCR(r)') = CRmax;
	
	% pbest
	pbest = 1 + floor(max(2, round(p * NP)) * rand(1, NP));
	
	% Archive
	XA = [X, A];
	SPA = [SP, A];
	
	% Index selection
	r1 = zeros(1, NP);
	r2 = zeros(1, NP);
	
	for i = 1 : NP		
		% Generate r1
		r1(i) = floor(1 + NP * rand);
		while i == r1(i)
			r1(i) = floor(1 + NP * rand);
		end
		
		% Generate r2
		r2(i) = floor(1 + (NP + nA) * rand);
		while i == r1(i) || r1(i) == r2(i)
			r2(i) = floor(1 + (NP + nA) * rand);
		end
	end
	
	% Mutation
	for i = 1 : NP	
		if FC(i) <= Q
			V(:, i) = X(:, i) ...
				+ F(i) .* (X(:, pbest(i)) - X(:, i)) ...
				+ F(i) .* (X(:, r1(i)) - XA(:, r2(i)));
		else
			V(:, i) = SP(:, i) ...
				+ F(i) .* (SP(:, sortidxfSP(pbest(i))) - SP(:, i)) ...
				+ F(i) .* (SP(:, r1(i)) - SPA(:, r2(i)));
		end
	end
	
	C = (1 - cw) * C + cw * cov(X');
	[B, ~] = eig(C);
	XT = X;
	VT = V;
	UT = U;
	SPT = SP;
	
	for i = 1 : NP
		jrand = floor(1 + D * rand);
		if FC(i) <= Q
			if rand < ER(i)
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
				U(:, i) = B * UT(:, i);
			else
				% Binominal Crossover
				for j = 1 : D
					if rand < CR(i) || j == jrand
						U(j, i) = V(j, i);
					else
						U(j, i) = X(j, i);
					end
				end
			end
		else
			if rand < ER(i)
				% SPS+EIG framework
				XT(:, i) = B' * X(:, i);
				VT(:, i) = B' * V(:, i);
				SPT(:, i) = B' * SP(:, i);
				for j = 1 : D
					if rand < CR(i) || j == jrand
						UT(j, i) = VT(j, i);
					else
						UT(j, i) = SPT(j, i);
					end
				end
				U(:, i) = B * UT(:, i);
			else
				% SPS framework
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
	
	% Constraint violation
	if EpsilonMethod
		for i = 1 : NP
			clbu = lb - U(:, i);
			cubu = U(:, i) - ub;
			psai_u(i) = sum(clbu(clbu > 0)) + sum(cubu(cubu > 0));
			
			if ~isempty(nonlcon)
				[cu, cequ] = feval(nonlcon, U(:, i));
				countcon = countcon + 1;
				psai_u(i) = psai_u(i) + sum(cu(cu > 0)) + sum(cequ(cequ > 0));
			end
		end
	end
	
	% Selection
	FailedIteration = true;
	nS = 0;
	if ~EpsilonMethod
		for i = 1 : NP
			if fu(i) < fx(i)
				nS			= nS + 1;
				S_ER(nS)	= ER(i);
				S_F(nS)		= F(i);
				S_CR(nS)	= CR(i);
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
			else
				FC(i) = FC(i) + 1;
			end
		end
	else
		% Epsilon level comparisons
		for i = 1 : NP
			X_AND_U_IN_EPSILON = psai_u(i) < EpsilonValue && psai_x(i) < EpsilonValue;
			X_AND_U_EQUAL_EPSILON = psai_u(i) == psai_x(i);
			
			if ((X_AND_U_IN_EPSILON || X_AND_U_EQUAL_EPSILON) && fu(i) < fx(i)) || ...
					(~X_AND_U_IN_EPSILON && psai_u(i) < psai_x(i))

				nS			= nS + 1;
				S_ER(nS)	= ER(i);
				S_F(nS)		= F(i);
				S_CR(nS)	= CR(i);
				S_df(nS)	= abs(fu(i) - fx(i));
				X(:, i)		= U(:, i);
				fx(i)		= fu(i);
				psai_x(i)	= psai_u(i);
				
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
			else
				FC(i) = FC(i) + 1;				
			end
		end
	end
	
	% Update MER, MF, and MCR
	if nS > 0
		w = S_df(1 : nS) ./ sum(S_df(1 : nS));
		
		if all(S_ER(1 : nS) == 0)
			MER(iM) = 0;
		else
			MER(iM) = sum(w .* S_ER(1 : nS));
			MER(iM) = max(-1, MER(iM));
			MER(iM) = min(2, MER(iM));
		end
		
		S_CR(S_CR == CRmin) = 0;
		S_CR(S_CR == CRmax) = 1;
		MCR(iM) = sum(w .* S_CR(1 : nS));
		MCR(iM) = max(0.05, MCR(iM));
		MCR(iM) = min(0.95, MCR(iM));
		
		MF(iM) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));		
		iM = mod(iM, H) + 1;
	end
	
	% Update cw
	cw = (1 - counteval / maxfunevals) * cwinit;
	
	% Sort	
	if ~EpsilonMethod
		[fx, fidx] = sort(fx);
		X = X(:, fidx);
		FC = FC(fidx);
	else
		PsaiFx = [psai_x', fx'];
		[~, SortingIndex] = sortrows(PsaiFx);
		X = X(:, SortingIndex);
		fx = fx(SortingIndex);
		FC = FC(SortingIndex);
		psai_x = psai_x(SortingIndex);
	end	
	
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
		'countcon', countcon, ...
		'muMF', mean(MF), ...
		'muMCR', mean(MCR), ...
		'muMER', mean(MER), ...
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
final.MER		= MER;
final.MCR		= MCR;
final.MF		= MF;
final.psai		= psai_x;
final.iM		= iM;
final.FC		= FC;
final.C			= C;
final.SP		= SP;
final.fSP		= fSP;
final.iSP		= iSP;
final.counteval = counteval;
final.countiter = countiter;
final.countstagnation = countstagnation;
final.countcon	= countcon;

out = finishoutput(out, X, fx, counteval, countiter, ...
	'countcon', countcon, ...
	'final', final, ...
	'muMF', mean(MF), ...
	'muMCR', mean(MCR), ...
	'muMER', mean(MER), ...
	'muFC', mean(FC));

% out = finishoutput(out, X, fx, counteval, countiter, ...
% 	'countcon', countcon, ...
% 	'muMF', mean(MF), ...
% 	'muMCR', mean(MCR), ...
% 	'muMER', mean(MER), ...
% 	'muFC', mean(FC));
end
