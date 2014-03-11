function [xmin, fmin, out] = apderand1bin(fitfun, lb, ub, maxfunevals, options)
% DERAND1BIN DE/rand/1/bin with adaptive poplution size 
% DERAND1BIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DERAND1BIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 0;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = 0;
defaultOptions.TolX = 0;
defaultOptions.TolStagnationIteration = 20;
defaultOptions.initPopFac = 5;
defaultOptions.MovingAverage = 20;
defaultOptions.Ptarget = 1e-4;
defaultOptions.NP_MIN = 48;
defaultOptions.PopFac = 0.98;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];

defaultOptions.TolCon = 1e-6;
defaultOptions.nonlcon = [];
defaultOptions.initial.cm = []; % Constraint violation measure
defaultOptions.initial.nc = []; % Number of violated constraints

options = setdefoptions(options, defaultOptions);
CR = options.CR;
F = options.F;
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
NP_MIN = options.NP_MIN;
initPopFac = options.initPopFac;
PopFac = options.PopFac;

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	f = options.initial.f;
	cm = options.initial.cm;
	nc = options.initial.nc;
else
	X = [];
	f = [];
	cm = [];
	nc = [];
end

D = numel(lb);
if isempty(X)
	if D == 10
		T = 0.9674;
	elseif D == 30
		T = 0.9874;
	else
		T = 0.9885;
	end
	
	if D == 10 || D == 30 || D == 50
		NP = floor(maxfunevals / ...
			(log(Ptarget) - log(sqrt(mean(ub - lb)^2/12))) * log(T));
		NP = max(NP, NP_MIN);
	else
		NP = ceil(D * initPopFac);
	end
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
	'converg_rate');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = preparecontourdata(D, lb, ub, fitfun);
end

% Initialize population
if isempty(X)
	if NP < 1e1
		LHS = lhsdesign(NP, D, 'iteration', 10)';
	elseif NP < 1e2
		LHS = lhsdesign(NP, D, 'iteration', 2)';
	else
		LHS = rand(D, NP);
	end
	
	X = zeros(D, NP);
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* LHS(:, i);
	end
end

% Constraint violation measure
if isempty(cm) || isempty(nc)
	cm = zeros(1, NP);
	nc = zeros(1, NP);
	
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
	f = zeros(1, NP);
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
nf = f;
nf(isinf(nf)) = [];
nfmax = max(nf);
nfmin = min(nf);
ncm = cm;
ncmmax = max(ncm);
ncmmin = min(ncm);

for i = 1 : NP
	if nc(i) == 0
		pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
	else
		pf(i) = nc(i) + (ncm(i) - ncmmin) / (ncmmax - ncmmin + eps);
	end
end

[pf, pfidx] = sort(pf);
f = f(pfidx);
X = X(:, pfidx);
cm = cm(pfidx);
nc = nc(pfidx);

% Initialize variables
V = X;
U = X;
cm_u = cm;
nc_u = nc;

% Display
if isDisplayIter
	displayitermessages(X, U, f, countiter, XX, YY, ZZ);
end

% Convergence speed
xstd_new = std(X, 0, 2);
m = 0;
Sconvrate = zeros(1, M);
SconvrateCounter = 0;
converg_rate = 1;

% Record
out = updateoutput(out, X, f, counteval, ...
	'NP', NP, ...
	'converg_rate', converg_rate);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(f) <= ftarget;
	fitnessconvergence = all(nc == 0) && isConverged(f, TolFun) ...
		&& isConverged(cm, TolCon);
	solutionconvergence = isConverged(X, TolX);
	stagnation = countStagnation >= TolStagnationIteration;
	
	% Convergence conditions
	if outofmaxfunevals || reachftarget || fitnessconvergence || ...
			solutionconvergence || stagnation
		break;
	end
	
	% Mutation
	for i = 1 : NP
		r1 = floor(1 + NP * rand);
		r2 = floor(1 + NP * rand);
		r3 = r2;
		
		while r2 == r3
			r3 = floor(1 + NP * rand);
		end
		
		V(:, i) = X(:, r1) + F * (X(:, r2) - X(:, r3));
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		for j = 1 : D
			if rand < CR || j == jrand
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
		displayitermessages(X, U, f, countiter, XX, YY, ZZ, ...
			'm', m);
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
			FailedIteration = false;
		end
	end
	
	% Convergence speed
	xstd_prev = xstd_new;
	xstd_new = std(X, 0, 2);
	converg_rate = mean(xstd_new ./ (xstd_prev + eps));
	
	if SconvrateCounter < M
		Sconvrate(SconvrateCounter + 1) = converg_rate;		
		SconvrateCounter = SconvrateCounter + 1;
	else
		Sconvrate(1:end-1) = Sconvrate(2:end);
		Sconvrate(end) = converg_rate;
		convrate_avg = geomean(Sconvrate);
		m = (log(Ptarget) - log(mean(xstd_new))) / (log(convrate_avg) + eps);
		
		if m > 0 && m > (maxfunevals - counteval) / NP
			NP_prev = NP;
			NP = max(floor(PopFac * NP_prev), NP_MIN);
			randidx = randperm(NP_prev, NP);
			X = X(:, randidx);
			V = V(:, randidx);
			U = U(:, randidx);
			f = f(randidx);
			pf = pf(randidx);
			cm = cm(randidx);
			nc = nc(randidx);
		end
	end
	
	% Sort
	nf = f;
	nf(isinf(nf)) = [];
	nfmax = max(nf);
	nfmin = min(nf);
	ncm = cm;
	ncmmax = max(ncm);
	ncmmin = min(ncm);
	
	for i = 1 : NP
		if nc(i) == 0
			pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
		else
			pf(i) = nc(i) + (ncm(i) - ncmmin) / (ncmmax - ncmmin + eps);
		end
	end
	
	[pf, pfidx] = sort(pf);
	f = f(pfidx);
	X = X(:, pfidx);
	cm = cm(pfidx);
	nc = nc(pfidx);
			
	% Record
	out = updateoutput(out, X, f, counteval, ...
		'NP', NP, ...
		'converg_rate', converg_rate);
	
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

final.cm = cm;
final.nc = nc;

out = finishoutput(out, X, f, counteval, ...
	'final', final, ...
	'countcon', countcon, ...
	'NP', NP, ...
	'converg_rate', converg_rate);
end