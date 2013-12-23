function [xbest1, xbest2, fbest, out] = mmjade_o(fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, options1, options2)
% MMJADE_O Sequential Min-max JADE
% MMJADE_O(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MMJADE_O(..., options1) minimizes the function with the given
% options Options1 for the 1st layer.
% MMJADE_O(..., options1, options2) minimize the function with the
% given options Options2 for the 2nd layer.

% Check number of input arguments
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
defaultOptions1.delta_F = 0.1;
defaultOptions1.delta_CR = 0.1;
defaultOptions1.p = 0.05;
defaultOptions1.w = 0.1;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolX = 0;
defaultOptions1.TolFun = 0;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.InnerSolver = 'jadebin';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.A = [];
defaultOptions1.initial.mu_F = [];
defaultOptions1.initial.mu_CR = [];
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
defaultOptions2.CR = 0.5;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.TolFun = 0;
defaultOptions2.TolX = 0;
options2 = setdefoptions(options2, defaultOptions2);

% Initialize algorithmic variables
dimensionFactor = max(1, options1.dimensionFactor);
delta_F = options1.delta_F;
delta_CR = options1.delta_CR;
p = options1.p;
w = options1.w;
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
	options1.initial = setdefoptions(options1.initial, defaultOptions1.initial);
	X = options1.initial.X;
	f = options1.initial.f;
	A = options1.initial.A;
	mu_F = options1.initial.mu_F;
	mu_CR = options1.initial.mu_CR;
	cm = options1.initial.cm;
	nc = options1.initial.nc;
	innerState = options1.initial.innerState;
else
	X = [];
	f = [];
	A = [];
	mu_F = [];
	mu_CR = [];
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

% Initialize archive
if isempty(A)
	A = zeros(D1, 2 * NP1);
	A(:, 1 : NP1) = X;
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
V = X;
U = X;
PX2 = zeros(D2, NP2, NP1);
pbest_size = p * NP1;
fu = zeros(1, NP1);
innerOutX = cell(1, NP1);
innerOutU = cell(1, NP1);
cm_u = zeros(1, NP1);
nc_u = zeros(1, NP1);

out = initoutput(RecordPoint, D1, NP1, maxfunevals, ...
	'innerFstd', ...
	'innerMeanXstd', ...
	'successRate', ...
	'X_Converged_FEs', ...
	'U_Converged_FEs', ...
	'mu_F', ...
	'mu_CR');

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

% mu_F
if isempty(mu_F)
	mu_F = options1.F;
end

% mu_CR
if isempty(mu_CR)
	mu_CR = options1.CR;
end

% Display
if isDisplayIter
	if all(isinf(f))
		dispconitermsg([X; innerXbest], [U; innerUbest], ...
			cm, countiter, ...
			XX, YY, ZZ, CC, 'counteval', counteval, ...
			'successRate', successRate, ...
			'mu_F', mu_F, ...
			'mu_CR', mu_CR);
	else
		dispconitermsg([X; innerXbest], [U; innerUbest], ...
			f(~isinf(f)), countiter, ...
			XX, YY, ZZ, CC, 'counteval', counteval, ...
			'successRate', successRate, ...
			'mu_F', mu_F, ...
			'mu_CR', mu_CR);
	end
	
	display_inner_info(innerState);
end

% Record minimal function values
out = updateoutput(out, X, f, counteval, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'mu_F', mu_F, ...
	'mu_CR', mu_CR);

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
	
	% Scaling factor and crossover rate
	S_CR = zeros(1, NP1);
	CR = mu_CR + delta_CR * randn(1, NP1);
	CR(CR > 1) = 1 - eps;
	CR(CR < 0) = eps;
	S_F = zeros(1, NP1);
	F = cauchyrnd(mu_F, delta_F, NP1, 1);
	F(F > 1) = 1 - eps - 0.01 * rand;
	
	for retry = 1 : 3
		if all(F > 0)
			break;
		end
		
		F(F <= 0) = cauchyrnd(mu_F, delta_F, sum(F <= 0), 1);
		F(F > 1) = 1 - eps - 0.01 * rand;
	end
	
	F(F <= 0) = 0.01 * mu_F * (1 + rand);
	
	Succ_Counter = 0;
	XA = [X, A];
	
	% Mutation
	for i = 1 : NP1
		% Try generating V within bounds
		for retry_within_bounds = 1 : NP1
			
			% Generate pbest_idx
			for retry = 1 : 3
				pbest_idx = max(1, ceil(rand * pbest_size));
				if ~all(X(:, pbest_idx) == X(:, i))
					break;
				end
			end
			
			% Generate r1
			for retry = 1 : NP1
				r1 = floor(1 + NP1 * rand);
				if i ~= r1
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : NP1 * NP1
				r2 = floor(1 + 2 * NP1 * rand);
				if ~(all(X(:, i) == XA(:, r2)) || all(X(:, r1) == XA(:, r2)))
					break;
				end
			end
			
			% Generate Vi
			V(:, i) = X(:, i) + F(i) .* ...
				(X(:, pbest_idx) - X(:, i) + X(:, r1) - XA(:, r2));
			
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
			if rand < CR(i) || j == jrand
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
			
% 			% Copy from innerXbest
% 			n_migration = ceil(migrateFactor * NP2);
% 			beginIndex = NP2 - n_migration;
% 			for j = 1 : n_migration
% 				r = floor(NP1 * rand + 1);
% 				PX2(:, beginIndex + j, i) = innerXbest(:, r) .* (1 + zeta * randn(D2, 1));
% 			end
		end
	else		
		for i = 1 : NP1
			for j = 1 : NP2
				PX2(:, j, i) = lb2 + (ub2 - lb2) .* rand(D2, 1);
			end
		end
	end
	
	PU2 = PX2;
	
% 	if ~isempty(nonlcon)
% 		parfor i = 1 : NP1
% 			% Compute XC
% 			innerFitfunXCi = @(y) max(feval(nonlcon, X(:, i), y)).^2;
% 			
% 			innerOptionsXCi = options2;
% 			innerOptionsXCi.initial = [];
% 			
% 			[innerXCbest(:, i), ~, innerOutXC{i}] = ...
% 				feval(innerSolver, innerFitfunXCi, ...
% 				lb2, ub2, ...
% 				innerMaxfunevalsX, innerOptionsXCi); 
% 						
% 			% Compute UC
% 			innerFitfunUCi = @(y) max(feval(nonlcon, U(:, i), y)).^2;
% 			
% 			innerOptionsUCi = options2;
% 			innerOptionsUCi.initial = [];
% 			
% 			[innerUCbest(:, i), ~, innerOutUC{i}] = ...
% 				feval(innerSolver, innerFitfunUCi, ...
% 				lb2, ub2, ...
% 				innerMaxfunevalsX, innerOptionsUCi); 
% 		end
% 		
% 		for i = 1 : NP1
% 			PX2(:, end, i) = innerXCbest(:, i);
% 			PU2(:, end, i) = innerUCbest(:, i);
% 			
% 			countcon = countcon + innerOutXC{i}.fes(end);
% 			countcon = countcon + innerOutUC{i}.fes(end);
% 		end
% 	end
	
	% Selection
	innerMaxfunevalsX = innerMaxIter * NP2;
	parfor i = 1 : NP1
% 		% Compute fxi, f(i)
% 		innerFitfunXi = @(X2) -feval(fitfun, X(:, i), X2);
% 		optionsX2i = options2;
% 		optionsX2i.initial = [];
% 		optionsX2i.initial.X = PX2(:, :, i);
% 		
% 		if ~isempty(nonlcon)
% 			optionsX2i.nonlcon = @(X2) feval(nonlcon, X(:, i), X2);
% 		end
% 		
% 		[innerXbest(:, i), innerFbest, innerOutX{i}] = ...
% 			feval(innerSolver, innerFitfunXi, ...
% 			lb2, ub2, ...
% 			innerMaxfunevalsX, optionsX2i);
		
		X_Converged_FEs(i) = 0;
% 		f(i) = -innerFbest;
		
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
			A(:, NP1 + Succ_Counter + 1) = U(:, i);
			innerXbest(:, i) = innerUbest(:, i);
			innerState{i} = innerOutU{i}.final;
			successRate = successRate + 1 / NP1;
			S_F(Succ_Counter + 1) = F(i);
			S_CR(Succ_Counter + 1) = CR(i);
			Succ_Counter = Succ_Counter + 1;
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
				'successRate', successRate, ...
				'mu_F', mu_F, ...
				'mu_CR', mu_CR);
		else
			dispconitermsg([X; innerXbest], [U; innerUbest], ...
				f(~isinf(f)), countiter, ...
				XX, YY, ZZ, CC, 'counteval', counteval, ...
				'successRate', successRate, ...
				'mu_F', mu_F, ...
				'mu_CR', mu_CR);
		end
		
		display_inner_info(innerState);
	end
	
	% Update archive
	rand_idx = randperm(NP1 + Succ_Counter);
	A(:, 1 : NP1) = A(:, rand_idx(1 : NP1));
	
	% Update CR and F
	if Succ_Counter > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : Succ_Counter));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : Succ_Counter).^2) / sum(S_F(1 : Succ_Counter));
	else
		mu_F = (1-w) * mu_F;
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
	out = updateoutput(out, X, f, counteval, ...
		'innerFstd', computeInnerFstd(innerState), ...
		'innerMeanXstd', computeInnerMeanXstd(innerState), ...
		'successRate', successRate, ...
		'X_Converged_FEs', mean(X_Converged_FEs), ...
		'U_Converged_FEs', mean(U_Converged_FEs), ...
		'mu_F', mu_F, ...
		'mu_CR', mu_CR);
	
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

final.A = A;
final.mu_F = mu_F;
final.mu_CR = mu_CR;
final.cm = cm;
final.nc = nc;
final.innerState = innerState;

out = finishoutput(out, X, f, counteval, 'final', final, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'mu_F', mu_F, ...
	'mu_CR', mu_CR, ...
	'countcon', countcon);
end