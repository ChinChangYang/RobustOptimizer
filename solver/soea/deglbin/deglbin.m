function [xmin, fmin, out] = deglbin(fitfun, lb, ub, maxfunevals, options)
% DEGLBIN Differential Evolution with global and local neighborhoods
% DEGLBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DEGLBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 10;
defaultOptions.CR = 0.5;
defaultOptions.NeighborhoodRatio = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.TolStagnationIteration = 20;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.w = [];

defaultOptions.TolCon = 1e-6;
defaultOptions.nonlcon = [];
defaultOptions.initial.cm = []; % Constraint violation measure
defaultOptions.initial.nc = []; % Number of violated constraints

options = setdefoptions(options, defaultOptions);
dimensionFactor = max(1, options.dimensionFactor);
CR = options.CR;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
TolStagnationIteration = options.TolStagnationIteration;
TolCon = options.TolCon;
nonlcon = options.nonlcon;

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	f = options.initial.f;
	w = options.initial.w;
	cm = options.initial.cm;
	nc = options.initial.nc;
else
	X = [];
	f = [];
	w = [];
	cm = [];
	nc = [];
end

D = numel(lb);

if isempty(X)
	NP = ceil(dimensionFactor * D);
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countcon = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals);

if isempty(w)	
	w = 0.05 + 0.9 * rand(1, NP);
end

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
	
	% dirty magic	
	X(:, end-1) = lb + abs(1e-7 * randn);
	X(:, end) = ub - abs(1e-7 * randn);
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
w = w(pfidx);

% Initialize variables
k = ceil(0.5 * (options.NeighborhoodRatio * NP));
wc = w;
V = X;
U = X;
cm_u = cm;
nc_u = nc;

% Display
if isDisplayIter
	displayitermessages(X, U, f, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, f, counteval);

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
	% Global best
	g_best = 1;
	
	for i = 1 : NP
		% Generate random mutant factor F, and parameters, alpha and beta.
		F = abs(0.5 * log(rand));
		alpha = F;
		beta = F;
		
		% Neiborhoods index
		n_index = (i-k) : (i+k);
		lessthanone = n_index < 1;
		n_index(lessthanone) = n_index(lessthanone) + NP;
		greaterthanNP = n_index > NP;
		n_index(greaterthanNP) = n_index(greaterthanNP) - NP;
		
		% Neiborhood solutions and fitness
		Xn = X(:, n_index);
		fn = f(n_index);
		
		% Best neiborhood
		[~, n_besti] = min(fn);
		Xn_besti = Xn(:, n_besti);
		
		% Random neiborhood index
		n_index(n_index == i) = [];
		Xn = X(:, n_index);
		p = ceil(rand * numel(n_index));
		q = ceil(rand * numel(n_index));
		
		while p == q
			q = ceil(rand * numel(n_index));
		end
		
		% Random neiborhood solutions
		Xp = Xn(:, p);
		Xq = Xn(:, q);
		
		% Local donor vector
		Li = X(:, i) + alpha * (Xn_besti - X(:, i)) + ...
			beta * (Xp - Xq);
		
		% Global donor vector
		r1 = floor(1 + NP * rand);
		
		while i == r1
			r1 = floor(1 + NP * rand);
		end
		
		r2 = floor(1 + NP * rand);
		
		while i == r2 || r1 == r2
			r2 = floor(1 + NP * rand);
		end
		
		gi = X(:, i) + alpha * (X(:, g_best) - X(:, i)) + ...
			beta * (X(:, r1) - X(:, r2));
		
		% Self-adaptive weight factor
		wc(i) = w(i) + F * (w(g_best) - w(i)) + ...
			F * (w(r1) - w(r2));
		
		if wc(i) < 0.05
			wc(i) = 0.05;
		elseif wc(i) > 0.95
			wc(i) = 0.95;
		end
		
		V(:, i) = wc(i) * gi + (1 - wc(i)) * Li;
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
	
	% Display
	if isDisplayIter
		displayitermessages(X, U, f, countiter, XX, YY, ZZ);
	end
	
	% Repair
	for i = 1 : NP
		for j = 1 : D
			if U(j, i) < lb(j)
				U(j, i) = lb(j);
			elseif U(j, i) > ub(j)
				U(j, i) = ub(j);
			end
		end
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
			w(i) = wc(i);
			FailedIteration = false;
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
	w = w(pfidx);
	
	% Record
	out = updateoutput(out, X, f, counteval);
	
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

% Final state of this algorithm
final.cm = cm;
final.nc = nc;
final.w = w;

out = finishoutput(out, X, f, counteval, ...
	'final', final, ...
	'countcon', countcon);
end