function [xbest1, xbest2, fbest, out] = minmaxdegl(fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXDEGL A differential evolution approach for solving constrained
% min-max optimization problems in Expert Systems with Applications 39
% (2013) 13440-13450 by Gilberto A.S. Segundo, Renato A. Krohling, and
% Rodrigo C. Cosme.
% MINMAXDEGL(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MINMAXDEGL(..., options1) minimizes the function with the given options
% Options1 for the 1st layer.
% MINMAXDEGL(..., options1, options2) minimize the function with the given
% options Options2 for the 2nd layer.

% Input arguments
if nargin <= 6
	options1 = [];
end

if nargin <= 7
	options2 = [];
end

D1 = numel(lb1);
D2 = numel(lb2);

% Default options for Layer 1
defaultOptions1.dimensionFactor = 30;
defaultOptions1.CR = 0.9;
defaultOptions1.NeighborhoodRatio = 0.1;
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolX = 0;
defaultOptions1.TolFun = 0;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.innerMaxIter = 200;
defaultOptions1.InnerSolver = 'deglbin';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.cm = [];
defaultOptions1.initial.nc = [];
defaultOptions1.initial.innerState = [];

defaultOptions1.TolCon = 1e-6;
defaultOptions1.nonlcon = [];
options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2.dimensionFactor = 30;
defaultOptions2.CR = 0.9;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.TolFun = 0;
defaultOptions2.TolX = 0;
options2 = setdefoptions(options2, defaultOptions2);

% Initialize algorithmic variables
dimensionFactor = max(1, options1.dimensionFactor);
CR = options1.CR;
isDisplayIter = strcmp(options1.Display, 'iter');
RecordPoint = max(0, floor(options1.RecordPoint));
TolFun = options1.TolFun;
TolX = options1.TolX;
TolStagnationIteration = options1.TolStagnationIteration;
innerSolver = options1.InnerSolver;
TolCon = options1.TolCon;
nonlcon = options1.nonlcon;
innerMaxIter = options1.innerMaxIter;

if ~isempty(options1.initial)
	options.initial = setdefoptions(options1.initial, defaultOptions1.initial);
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

innerState = options1.initial.innerState;
existInnerState = ~isempty(innerState);

NP1 = ceil(dimensionFactor * D1);
NP2 = ceil(options2.dimensionFactor * D2);

nBest = round(0.2 * NP1);

% Initialize contour data
if isDisplayIter
	contourOptions.nonlcon = nonlcon;
	[XX, YY, ZZ, CC] = ...
		cminmaxcontourdata(D1, lb1, ub1, lb2, ub2, fitfun, contourOptions);
end

% Initialize population
if isempty(X)
	if NP1 < 1e1
		LHS = lhsdesign(NP1, D1, 'iteration', 10)';
	elseif NP1 < 1e2
		LHS = lhsdesign(NP1, D1, 'iteration', 2)';
	else
		LHS = rand(D1, NP1);
	end
	
	X = zeros(D1, NP1);
	for i = 1 : NP1
		X(:, i) = lb1 + (ub1 - lb1) .* LHS(:, i);
	end
end

% Initialize inner states
if isempty(innerState)
	innerState = cell(1, NP1);
end

% Initialize variables
counteval = 0;
countcon = 0;
countiter = 1;
countStagnation = 0;
successRate = 0;
U_Converged_FEs = zeros(1, NP1);
innerXbest = zeros(D2, NP1);
innerUbest = innerXbest;
V = X;
U = X;
V2 = zeros(D2, NP2, NP1);
fu = zeros(1, NP1);
innerOutX = cell(1, NP1);
innerOutU = cell(1, NP1);
cm_u = zeros(1, NP1);
nc_u = zeros(1, NP1);
k = ceil(0.5 * (options1.NeighborhoodRatio * NP1));
w = 0.05 + 0.9 * rand(1, NP1);
wc = w;
out = initoutput(RecordPoint, D1, NP1, maxfunevals, ...
	'innerFstd', ...
	'innerMeanXstd', ...
	'successRate', ...
	'U_Converged_FEs');

% Evaluation
if isempty(f)
	f = zeros(1, NP1);
	innerMaxfunevalsX = innerMaxIter * NP2;
	
	for i = 1 : NP1
		innerFitfun = @(y) -feval(fitfun, X(:, i), y);
		optionsX2i = options2;
		
		if ~isempty(nonlcon)
			innerNonlcon = @(y) feval(nonlcon, X(:, i), y);
			optionsX2i.nonlcon = innerNonlcon;
		end
		
		if existInnerState
			optionsX2i.initial = innerState{i};
		end
		
		[innerXbest(:, i), innerFbest, innerOutX{i}] = ...
			feval(innerSolver, innerFitfun, ...
			lb2, ub2, innerMaxfunevalsX, optionsX2i);
		
		f(i) = -innerFbest;
		innerState{i} = innerOutX{i}.final;
	end
	
	for i = 1 : NP1		
		counteval = counteval + innerOutX{i}.fes(end);
		countcon = countcon + innerOutX{i}.countcon;
	end
end

% Constraint violation measure
if isempty(cm) || isempty(nc)
	cm = zeros(1, NP1);
	nc = zeros(1, NP1);
	
	for i = 1 : NP1
		clb = lb1 - X(:, i);
		cub = X(:, i) - ub1;
		cm(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = sum(clb > 0) + sum(cub > 0);
	end
	
	for i = 1 : NP1
		clb = lb2 - innerXbest(:, i);
		cub = innerXbest(:, i) - ub2;
		cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
	end
	
	if ~isempty(nonlcon)
		for i = 1 : NP1
			[cx, ceqx] = feval(nonlcon, X(:, i), innerXbest(:, i));
			countcon = countcon + 1;
			cm(i) = cm(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
			nc(i) = nc(i) + sum(cx > 0) + sum(ceqx > 0);
		end
	end
end

% Sort
pf = zeros(1, NP1);
nf = f;
nf(isinf(nf)) = [];
nfmax = max(nf);
nfmin = min(nf);
cmmax = max(cm);
cmmin = min(cm);

for i = 1 : NP1
	if nc(i) == 0
		pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
	else
		pf(i) = nc(i) + (cm(i) - cmmin) / (cmmax - cmmin + eps);
	end
end

[pf, pfidx] = sort(pf);
f = f(pfidx);
X = X(:, pfidx);
innerXbest = innerXbest(:, pfidx);
innerState = innerState(pfidx);
cm = cm(pfidx);
nc = nc(pfidx);

% Display
if isDisplayIter
	if all(isinf(f))
		dispconitermsg([X; innerXbest], [U; innerUbest], ...
			cm, countiter, ...
			XX, YY, ZZ, CC, 'counteval', counteval, ...
			'successRate', successRate);
	else
		dispconitermsg([X; innerXbest], [U; innerUbest], ...
			f(~isinf(f)), countiter, ...
			XX, YY, ZZ, CC, 'counteval', counteval, ...
			'successRate', successRate);
	end
	
	display_inner_info(innerState);
	
	retry_print = true;
	while retry_print
		try
			filename = sprintf('minmaxdegl_%s_%d.eps', fitfun, countiter);
			print(filename, '-depsc');
			retry_print = false;
		catch ME
			if strcmp(ME.identifier, 'MATLAB:fileio:cantOpenFileNoPermission')
				retry_print = true;
			else
				rethrow(ME);
			end
		end
	end
end

% Record
out = updateoutput(out, X, f, counteval + countcon, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'U_Converged_FEs', mean(U_Converged_FEs));

countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval + countcon >= maxfunevals;
	fitnessconvergence = isConverged(f, TolFun) && isConverged(cm, TolCon);
	solutionconvergence = isConverged(X, TolX);
	stagnation = countStagnation >= TolStagnationIteration;
	
	% Convergence conditions
	if outofmaxfunevals || fitnessconvergence || solutionconvergence ...
			|| stagnation
		break;
	end
	
	% Global best solution
	[~, g_best] = min(pf);
	
	for i = 1 : NP1
		% Generate random mutant factor F, and parameters, alpha and beta.
		F = abs(0.5 * log(rand));
		alpha = F;
		beta = F;
		
		% Neiborhoods index
		n_index = (i-k) : (i+k);
		lessthanone = n_index < 1;
		n_index(lessthanone) = n_index(lessthanone) + NP1;
		greaterthanNP = n_index > NP1;
		n_index(greaterthanNP) = n_index(greaterthanNP) - NP1;
		
		% Neiborhood solutions and fitness
		Xn = X(:, n_index);
		pfn = pf(n_index);
		
		% Best neiborhood
		[~, n_besti] = min(pfn);
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
		r1 = floor(1 + NP1 * rand);
		
		while i == r1
			r1 = floor(1 + NP1 * rand);
		end
		
		r2 = floor(1 + NP1 * rand);
		
		while i == r2 || r1 == r2
			r2 = floor(1 + NP1 * rand);
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
		U(:, i) = V(:, i);
		
		% Crossover
		jrand = floor(1 + rand * D1);
		for j = 1 : D1
			if rand >= CR && j ~= jrand
				U(j, i) = X(j, i);
			end
		end
		
		% Calculate the distance
		d = zeros(1, NP1);
		for j = 1 : NP1
			d(j) = norm(U(:, i) - X(:, j));
		end
		
		% Sort distances
		[~, sort_d_index] = sort(d);
		bestPopC = innerXbest(:, sort_d_index(1 : nBest));
		
		for j = 1 : NP2
			V2(:, j, i) = lb2 + (ub2 - lb2) .* rand(D2, 1);
		end
		
		V2(:, 1 : nBest, i) = bestPopC;
		
		% dirty magic
		V2(:, end - 1, i) = lb2 + abs(1e-7 * rand);
		V2(:, end, i) = ub2 - abs(1e-7 * rand);
	end
	
	% Selection
	innerMaxfunevalsX = innerMaxIter * NP2;
	parfor i = 1 : NP1		
		% Compute fui
		fitfunU2i = @(U2) -feval(fitfun, U(:, i), U2);
		optionsU2i = options2;
		optionsU2i.initial = [];
		optionsU2i.initial.X = V2(:, :, i);
		optionsU2i.initial.f = [];
		optionsU2i.initial.w = [];
		optionsU2i.initial.cm = [];
		optionsU2i.initial.nc = [];
		
		if ~isempty(nonlcon)
			optionsU2i.nonlcon = @(U2) feval(nonlcon, U(:, i), U2);
		end
		
		[innerUbest(:, i), innerFbest, innerOutU{i}] = ...
			feval(innerSolver, fitfunU2i, ...
			lb2, ub2, ...
			innerMaxfunevalsX, optionsU2i);
		
		U_Converged_FEs(i) = innerOutU{i}.fes(end);		
		fu(i) = -innerFbest;
	end
	
	for i = 1 : NP1
		counteval = counteval + innerOutU{i}.fes(end);
		countcon = countcon + innerOutU{i}.countcon;
	end
	
	% Constraint violation measure
	for i = 1 : NP1
		clb = lb1 - X(:, i);
		cub = X(:, i) - ub1;
		cm(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = sum(clb > 0) + sum(cub > 0);
		
		clb = lb1 - U(:, i);
		cub = U(:, i) - ub1;
		cm_u(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = sum(clb > 0) + sum(cub > 0);
	end
	
	for i = 1 : NP1
		clb = lb2 - innerXbest(:, i);
		cub = innerXbest(:, i) - ub2;
		cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
				
		clb = lb2 - innerUbest(:, i);
		cub = innerUbest(:, i) - ub2;
		cm_u(i) = cm_u(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = nc_u(i) + sum(clb > 0) + sum(cub > 0);
	end
	
	if ~isempty(nonlcon)
		for i = 1 : NP1
			[cx, ceqx] = feval(nonlcon, X(:, i), innerXbest(:, i));
			cm(i) = cm(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
			nc(i) = nc(i) + sum(cx > 0) + sum(ceqx > 0);
			
			[cu, cequ] = feval(nonlcon, U(:, i), innerUbest(:, i));
			cm_u(i) = cm_u(i) + sum(cu(cu > 0)) + sum(cequ(cequ > 0));
			nc_u(i) = nc_u(i) + sum(cu > 0) + sum(cequ > 0);
		end
	end
	
	% Replacement
	successRate = 0;
	FailedIteration = true;
	for i = 1 : NP1
		if nc(i) == 0 && nc_u(i) == 0
			if fu(i) < f(i)
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
			f(i) = fu(i);
			X(:, i) = U(:, i);
			w(i) = wc(i);
			innerXbest(:, i) = innerUbest(:, i);
			innerState{i} = innerOutU{i}.final;
			successRate = successRate + 1 / NP1;
			FailedIteration = false;
		end
	end
	
	% Display
	if isDisplayIter
		if all(isinf(f))
			dispconitermsg([X; innerXbest], [U; innerUbest], ...
				cm, countiter, ...
				XX, YY, ZZ, CC, 'counteval', counteval, ...
				'successRate', successRate);
		else
			dispconitermsg([X; innerXbest], [U; innerUbest], ...
				f(~isinf(f)), countiter, ...
				XX, YY, ZZ, CC, 'counteval', counteval, ...
				'successRate', successRate);
		end
		
		display_inner_info(innerState);
		
		retry_print = true;
		while retry_print
			try
				filename = sprintf('minmaxdegl_%s_%d.eps', fitfun, countiter);
				print(filename, '-depsc');
				retry_print = false;
			catch ME
				if strcmp(ME.identifier, 'MATLAB:fileio:cantOpenFileNoPermission')
					retry_print = true;
				else
					rethrow(ME);
				end
			end
		end
	end
	
	% Sort
	nf = f;
	nf(isinf(nf)) = [];
	nfmax = max(nf);
	nfmin = min(nf);
	cmmax = max(cm);
	cmmin = min(cm);
	
	for i = 1 : NP1
		if nc(i) == 0
			pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
		else
			pf(i) = nc(i) + (cm(i) - cmmin) / (cmmax - cmmin + eps);
		end
	end
	
	[pf, pfidx] = sort(pf);
	cm = cm(pfidx);
	nc = nc(pfidx);
	f = f(pfidx);
	X = X(:, pfidx);
	w = w(pfidx);
	innerXbest = innerXbest(:, pfidx);
	innerState = innerState(pfidx);
	
	% Record
	out = updateoutput(out, X, f, counteval + countcon, ...
		'innerFstd', computeInnerFstd(innerState), ...
		'innerMeanXstd', computeInnerMeanXstd(innerState), ...
		'successRate', successRate, ...
		'U_Converged_FEs', mean(U_Converged_FEs));
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end
end

% The best individual
xbest1 = X(:, 1);
xbest2 = innerState{1}.X(:, 1);
fbest = f(1);

final.innerState = innerState;

out = finishoutput(out, X, f, counteval + countcon, ...
	'final', final, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'countcon', countcon);
end