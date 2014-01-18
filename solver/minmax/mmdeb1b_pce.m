function [xbest1, xbest2, fbest, out] = mmdeb1b_pce(fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, options1, options2)
% MMDEB1B_PCE Sequential Min-Max DE/best/1/bin with finding
% constraint boundary
% MMDEB1B_PCE(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MMDEB1B_PCE(..., options1) minimizes the function with the given
% options Options1 for the 1st layer.
% MMDEB1B_PCE(..., options1, options2) minimize the function with the
% given options Options2 for the 2nd layer.
if nargin <= 6
	options1 = [];
end

if nargin <= 7
	options2 = [];
end

% Determine dimension
D1 = numel(lb1);
D2 = numel(lb2);

% Default options for Layer 1
defaultOptions1.dimensionFactor = 10;
defaultOptions1.F = 0.9;
defaultOptions1.CR = 0.5;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolX = 0;
defaultOptions1.TolFun = 0;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.InnerSolver = 'debest1bin';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.cm = [];
defaultOptions1.initial.nc = [];
defaultOptions1.initial.innerState = [];

defaultOptions1.TolCon = 1e-6;
defaultOptions1.nonlcon = [];
defaultOptions1.innerMaxIter = 200;
defaultOptions1.migrateFactor = 0.7;

options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2.dimensionFactor = 10;
defaultOptions2.F = 0.9;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.TolFun = 0;
defaultOptions2.TolX = 0;
options2 = setdefoptions(options2, defaultOptions2);

% Initialize algorithmic variables
dimensionFactor = max(1, options1.dimensionFactor);
F = options1.F;
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
migrateFactor = options1.migrateFactor;

if ~isempty(options1.initial)
	options1.initial = setdefoptions(options1.initial, defaultOptions1.initial);
	X = options1.initial.X;
	f = options1.initial.f;
	cm = options1.initial.cm;
	nc = options1.initial.nc;
	innerState = options1.initial.innerState;
else
	X = [];
	f = [];
	cm = [];
	nc = [];
	innerState = [];
end

existInnerState = ~isempty(innerState);

NP1 = ceil(dimensionFactor * D1);
NP2 = ceil(options2.dimensionFactor * D2);

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
X_Converged_FEs = zeros(1, NP1);
U_Converged_FEs = zeros(1, NP1);
innerXbest = zeros(D2, NP1);
innerUbest = innerXbest;
innerXCbest = innerXbest;
innerUCbest = innerXbest;
V = X;
U = X;
PX2 = zeros(D2, NP2, NP1);
fu = zeros(1, NP1);
innerOutX = cell(1, NP1);
innerOutU = cell(1, NP1);
innerOutXC = cell(1, NP1);
innerOutUC = cell(1, NP1);
cm_u = zeros(1, NP1);
nc_u = zeros(1, NP1);

out = initoutput(RecordPoint, D1, NP1, maxfunevals, ...
	'innerFstd', ...
	'innerMeanXstd', ...
	'successRate', ...
	'X_Converged_FEs', ...
	'U_Converged_FEs');

% Evaluation
if isempty(f)
	f = zeros(1, NP1);
	innerMaxfunevalsX = innerMaxIter * NP2;
	
	parfor i = 1 : NP1
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
end

% Record minimal function values
out = updateoutput(out, X, f, counteval + countcon, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
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
	
	% Mutation
	for i = 1 : NP1
		% Try generating V within bounds
		for retry_within_bounds = 1 : NP1
			
			% Generate r1
			r1 = floor(1 + NP1 * rand);
			
			% Generate r2
			for retry = 1 : 3
				r2 = floor(1 + NP1 * rand);
				if ~(all(X(:, r1) == X(:, r2)))
					break;
				end
			end
			
			% Generate Fi
			Fi = F + 0.01 * randn;
			
			% Generate Vi
			V(:, i) = X(:, 1) + Fi .* (X(:, r1) - X(:, r2));
			
			% Check boundary
			if all(V(:, i) >= lb1) && all(V(:, i) <= ub1)
				break;
			end
		end
	end
	
	% Crossover
	for i = 1 : NP1
		jrand = floor(1 + D1 * rand);
		
		for j = 1 : D1
			if rand < CR || j == jrand
				U(j, i) = V(j, i);
			else
				U(j, i) = X(j, i);
			end
		end
	end
	
	% Prediction
	anyEmptyInnerState = false;
	for i = 1 : NP1
		if isempty(innerState{i})
			anyEmptyInnerState = true;
			break;
		end
	end
	
	if ~anyEmptyInnerState
		for i = 1 : NP1
			% Copy from itselfs individuals
			for j = 1 : NP2
				PX2(:, j, i) = innerState{i}.X(:, j);
			end
			
			% Copy from innerXbest
			n_migration = ceil(migrateFactor * NP2);
			beginIndex = NP2 - n_migration;
			for j = 1 : n_migration
				r = floor(NP1 * rand + 1);
				PX2(:, beginIndex + j, i) = innerXbest(:, r);
			end
		end
	else
		for i = 1 : NP1
			for j = 1 : NP2
				PX2(:, j, i) = lb2 + (ub2 - lb2) .* rand(D2, 1);
			end
		end
	end
	
	PU2 = PX2;
	
	if ~isempty(nonlcon)
		parfor i = 1 : NP1
			% Compute XC
			innerFitfunXCi = @(y) max(feval(nonlcon, X(:, i), y)).^2;
			
			innerOptionsXCi = options2;
			innerOptionsXCi.initial = [];
			
			[innerXCbest(:, i), ~, innerOutXC{i}] = ...
				feval(innerSolver, innerFitfunXCi, ...
				lb2, ub2, ...
				innerMaxfunevalsX, innerOptionsXCi); 
			
			% Compute UC
			innerFitfunUCi = @(y) max(feval(nonlcon, U(:, i), y)).^2;
			
			innerOptionsUCi = options2;
			innerOptionsUCi.initial = [];
			
			[innerUCbest(:, i), ~, innerOutUC{i}] = ...
				feval(innerSolver, innerFitfunUCi, ...
				lb2, ub2, ...
				innerMaxfunevalsX, innerOptionsUCi); 
		end
		
		for i = 1 : NP1
			PX2(:, end, i) = innerXCbest(:, i);
			PU2(:, end, i) = innerUCbest(:, i);
			
			countcon = countcon + innerOutXC{i}.fes(end);
			countcon = countcon + innerOutUC{i}.fes(end);
		end
	end
	
	% Selection
	innerMaxfunevalsX = innerMaxIter * NP2;
	parfor i = 1 : NP1
		% Compute fxi, f(i)
		innerFitfunXi = @(X2) -feval(fitfun, X(:, i), X2);
		optionsX2i = options2;
		optionsX2i.initial = [];
		optionsX2i.initial.X = PX2(:, :, i);
		
		if ~isempty(nonlcon)
			optionsX2i.nonlcon = @(X2) feval(nonlcon, X(:, i), X2);
		end
		
		[innerXbest(:, i), innerFbest, innerOutX{i}] = ...
			feval(innerSolver, innerFitfunXi, ...
			lb2, ub2, ...
			innerMaxfunevalsX, optionsX2i);
		
		X_Converged_FEs(i) = innerOutX{i}.fes(end);
		f(i) = -innerFbest;
		
		% Compute fui
		fitfunU2i = @(U2) -feval(fitfun, U(:, i), U2);
		optionsU2i = options2;
		optionsU2i.initial = [];
		optionsU2i.initial.X = PU2(:, :, i);
		
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
		counteval = counteval + innerOutX{i}.fes(end);
		countcon = countcon + innerOutX{i}.countcon;
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
			countcon = countcon + 1;
			cm(i) = cm(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
			nc(i) = nc(i) + sum(cx > 0) + sum(ceqx > 0);
			
			[cu, cequ] = feval(nonlcon, U(:, i), innerUbest(:, i));
			countcon = countcon + 1;
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
			innerXbest(:, i) = innerUbest(:, i);
			innerState{i} = innerOutU{i}.final;
			successRate = successRate + 1 / NP1;
			FailedIteration = false;
		else
			innerState{i} = innerOutX{i}.final;
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
	innerXbest = innerXbest(:, pfidx);
	innerState = innerState(pfidx);
	
	% Record
	out = updateoutput(out, X, f, counteval + countcon, ...
		'innerFstd', computeInnerFstd(innerState), ...
		'innerMeanXstd', computeInnerMeanXstd(innerState), ...
		'successRate', successRate, ...
		'X_Converged_FEs', mean(X_Converged_FEs), ...
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

xbest1 = X(:, 1);
xbest2 = innerState{1}.X(:, 1);
fbest = f(1);

final.innerState = innerState;

out = finishoutput(out, X, f, counteval + countcon, 'final', final, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'countcon', countcon);
end