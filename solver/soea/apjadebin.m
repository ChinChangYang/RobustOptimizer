function [xmin, fmin, out] = apjadebin(fitfun, lb, ub, maxfunevals, options)
% JADEBIN JADE algorithm with adaptive population size
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
defaultOptions.p = 0.05;
defaultOptions.w = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = 0;
defaultOptions.TolX = 0;
defaultOptions.TolStagnationIteration = 20;
defaultOptions.MovingAverage = 10;
defaultOptions.Ptarget = 1e-6;
defaultOptions.PopFac = 0.98;
defaultOptions.NP_MAX = 4000;
defaultOptions.r = 0.1;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.A = [];
defaultOptions.initial.mu_CR = [];
defaultOptions.initial.mu_F = [];

defaultOptions.TolCon = 1e-6;
defaultOptions.nonlcon = [];
defaultOptions.initial.cm = []; % Constraint violation measure
defaultOptions.initial.nc = []; % Number of violated constraints

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
delta_CR = options.delta_CR;
delta_F = options.delta_F;
p = options.p;
w = options.w;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
TolStagnationIteration = options.TolStagnationIteration;
TolCon = options.TolCon;
nonlcon = options.nonlcon;
M = options.MovingAverage;
Ptarget = options.Ptarget;
NP_MIN = numel(lb);
NP_MAX = options.NP_MAX;
r = options.r;

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	f = options.initial.f;
	A = options.initial.A;
	mu_CR = options.initial.mu_CR;
	mu_F = options.initial.mu_F;
	cm = options.initial.cm;
	nc = options.initial.nc;
else
	X = [];
	f = [];
	A = [];
	mu_CR = [];
	mu_F = [];
	cm = [];
	nc = [];
end

D = numel(lb);
if isempty(X)
	NP = max(1, ceil(dimensionFactor * D));
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countcon = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'NP', ...
	'converg_rate', ...
	'geomean_Sconvrate', ...
	'm');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
end

% Initialize population
if isempty(X)	
	X = zeros(D, NP_MAX);
	for i = 1 : NP_MAX
		X(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
end

% Initialize archive
if isempty(A)
	A = X;
end

% Constraint violation measure
if isempty(cm) || isempty(nc)
	cm = nan(1, NP_MAX);
	nc = nan(1, NP_MAX);
	
	for i = 1 : NP
		clb = lb - X(:, i);
		cub = X(:, i) - ub;
		cm(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = sum(clb > 0) + sum(cub > 0);
	end
	
	if ~isempty(nonlcon)		
		for i = 1 : NP
			[c, ceq] = feval(nonlcon, X(:, i));
			countcon = countcon + 1;
			cm(i) = cm(i) + sum(c(c > 0)) + sum(ceq(ceq > 0));
			nc(i) = nc(i) + sum(c > 0) + sum(ceq > 0);
		end
	end
end

% Evaluation
if isempty(f)
	f = nan(1, NP_MAX);
	for i = 1 : NP
		if nc(i) > 0
			f(i) = inf;
		else
			f(i) = feval(fitfun, X(:, i));
			counteval = counteval + 1;
		end
	end
end

% Sort
pf = zeros(1, NP);
nf = f(1 : NP);
nf(isinf(nf)) = [];
nfmax = max(nf);
nfmin = min(nf);
ncm = cm(1 : NP);
ncmmax = max(ncm);
ncmmin = min(ncm);

for i = 1 : NP
	if nc(i) == 0
		pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
	else
		pf(i) = nc(i) + (ncm(i) - ncmmin) / (ncmmax - ncmmin + eps);
	end
end

[~, pfidx] = sort(pf);
f(1 : NP) = f(pfidx);
X(:, 1 : NP) = X(:, pfidx);
cm(1 : NP) = cm(pfidx);
nc(1 : NP) = nc(pfidx);

% mu_F
if isempty(mu_F)
	mu_F = options.F;
end

% mu_CR
if isempty(mu_CR)
	mu_CR = options.CR;
end

% Initialize variables
V = X(:, 1 : NP);
U = X(:, 1 : NP);
pbest_size = p * NP;
cm_u = cm(1 : NP);
nc_u = nc(1 : NP);

% Display
if isDisplayIter
	displayitermessages(...
		X(:, 1 : NP), U(:, 1 : NP), f(1 : NP), countiter, XX, YY, ZZ, 'mu_F', mu_F, 'mu_CR', mu_CR);
end

% Convergence speed
xstd_new = std(X(:, 1 : NP), 0, 2);
xstd_rec = xstd_new;
Sconvrate = ones(1, M);
SconvrateCounter = 0;
m = 0;

% Record
out = updateoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
	'NP', NP, ...
	'geomean_Sconvrate', geomean(Sconvrate), ...
	'm', m);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(f(1 : NP)) <= ftarget;
	fitnessconvergence = all(nc(1 : NP) == 0) && isConverged(f(1 : NP), TolFun) ...
		&& isConverged(cm(1 : NP), TolCon);
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
	
	% Mutation
	for i = 1 : NP		
		for checkRetry = 1 : 3			
			% Generate pbest_idx
			for retry = 1 : 3
				pbest_idx = max(1, ceil(rand * pbest_size));
				if ~all(X(:, pbest_idx) == X(:, i))
					break;
				end
			end
			
			% Generate r1
			for retry = 1 : 3
				r1 = floor(1 + NP * rand);
				if i ~= r1
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : 3
				r2 = floor(1 + 2 * NP * rand);
				if ~(all(X(:, i) == XA(:, r2)) || all(X(:, r1) == XA(:, r2)))
					break;
				end
			end
							
			V(:, i) = X(:, pbest_idx) + F(i) .* ...
				(X(:, i) - X(:, i) + X(:, r1) - XA(:, r2));
% 			V(:, i) = X(:, i) + F(i) .* ...
% 				(X(:, pbest_idx) - X(:, i) + X(:, r1) - XA(:, r2));
			
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
	
	% Constraint reflection
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
	
	% Constraint violation measure		
	for i = 1 : NP
		clb = lb - U(:, i);
		cub = U(:, i) - ub;
		cm_u(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = sum(clb > 0) + sum(cub > 0);
	end
	
	if ~isempty(nonlcon)		
		for i = 1 : NP
			[c, ceq] = feval(nonlcon, U(:, i));
			countcon = countcon + 1;
			cm_u(i) = cm_u(i) + sum(c(c > 0)) + sum(ceq(ceq > 0));
			nc_u(i) = nc_u(i) + sum(c > 0) + sum(ceq > 0);
		end
	end
	
	% Selection
	FailedIteration = true;
	for i = 1 : NP
		fui = inf;
		
		if nc(i) == 0 && nc_u(i) == 0
			fui = feval(fitfun, U(:, i));
			counteval = counteval + 1;
			
			if fui < f(i)
				u_selected = true;
			else
				u_selected = false;
			end
		elseif nc(i) > nc_u(i)
			u_selected = true;
		elseif nc(i) < nc_u(i)
			u_selected = false;
		else % nvc(i) == nvc_u(i) && nvc(i) ~= 0 && nvc_u(i) ~= 0
			if cm(i) > cm_u(i)
				u_selected = true;
			else
				u_selected = false;
			end
		end			
		
		if u_selected
			cm(i) = cm_u(i);
			nc(i) = nc_u(i);			
			f(i) = fui;
			X(:, i) = U(:, i);
			A(:, NP + A_Counter + 1) = U(:, i);
			S_CR(A_Counter + 1) = CR(i);
			S_F(A_Counter + 1) = F(i);
			A_Counter = A_Counter + 1;
			FailedIteration = false;
		end
	end
	
	% Update archive
	rand_idx = randperm(NP + A_Counter);
	A(:, 1 : NP) = A(:, rand_idx(1 : NP));
	
	% Update CR and F
	if A_Counter > 0
		mu_CR = (1 - w) * mu_CR + w * mean(S_CR(1 : A_Counter));
		mu_F = (1 - w) * mu_F + w * sum(S_F(1 : A_Counter).^2) / sum(S_F(1 : A_Counter));
	else
		mu_F = (1 - w) * mu_F;
	end
	
	% Sort	
	pf = zeros(1, NP);
	nf = f(1 : NP);
	nf(isinf(nf)) = [];
	nfmax = max(nf);
	nfmin = min(nf);
	ncm = cm(1 : NP);
	ncmmax = max(ncm);
	ncmmin = min(ncm);
	
	for i = 1 : NP
		if nc(i) == 0
			pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
		else
			pf(i) = nc(i) + (ncm(i) - ncmmin) / (ncmmax - ncmmin + eps);
		end
	end
	
	[~, pfidx] = sort(pf);
	f(1 : NP) = f(pfidx);
	X(:, 1 : NP) = X(:, pfidx);
	cm(1 : NP) = cm(pfidx);
	nc(1 : NP) = nc(pfidx);
	
	% Convergence speed
	xstd_prev = xstd_new;
	xstd_new = std(X(:, 1 : NP), 0, 2);	
	converg_rate = mean(xstd_new ./ (xstd_prev + eps));
	NP_prev = NP;
	
	if SconvrateCounter < 0
		SconvrateCounter = SconvrateCounter + 1;
	elseif SconvrateCounter < M
		Sconvrate(SconvrateCounter + 1) = converg_rate;		
		SconvrateCounter = SconvrateCounter + 1;
	else
		Sconvrate(1:end-1) = Sconvrate(2:end);
		Sconvrate(end) = converg_rate;
		convrate_avg = geomean(Sconvrate);
		m = (log(Ptarget) - log(mean(xstd_new))) / (log(convrate_avg) + eps);
		NP_exp = (maxfunevals - counteval) / m;
		
		if NP_exp < NP
			NP = floor((1 - r) * NP_prev + r * NP_exp);
			NP = max(NP_MIN, min(NP, NP_MAX));
			SconvrateCounter = -M;
		elseif mean(xstd_new) < mean(xstd_rec)
			xstd_rec = xstd_new;
			NP = floor((1 - r) * NP_prev + r * NP_exp);
			NP = max(NP_MIN, min(NP, NP_MAX));		
			SconvrateCounter = -M;
		end
	end
	
	% Check function values
	for i = NP_prev : NP
		if isnan(f(i))
			f(i) = feval(fitfun, X(:, i));
			
			clb = lb - X(:, i);
			cub = X(:, i) - ub;
			cm(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
			nc(i) = sum(clb > 0) + sum(cub > 0);
						
			if ~isempty(nonlcon)
				[c, ceq] = feval(nonlcon, X(:, i));
				countcon = countcon + 1;
				cm(i) = cm(i) + sum(c(c > 0)) + sum(ceq(ceq > 0));
				nc(i) = nc(i) + sum(c > 0) + sum(ceq > 0);
			end
		end
	end
	
	% Record
	out = updateoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
		'NP', NP, ...
		'geomean_Sconvrate', geomean(Sconvrate), ...
		'm', m);
	
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
final.cm = cm(1 : NP);
final.nc = nc(1 : NP);

out = finishoutput(out, X(:, 1 : NP), f(1 : NP), counteval, ...
	'final', final, ...
	'countcon', countcon, ...
	'NP', NP, ...
	'geomean_Sconvrate', geomean(Sconvrate), ...
	'm', m);
end
