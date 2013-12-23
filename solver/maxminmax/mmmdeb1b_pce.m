function [xbest1, xbest2, xbest3, fbest, out] = mmmdeb1b_c(...
	fitfun, ...		% objective function f(x,y,z)
	maxfunevals, ...	% maximal function evaluations
	lb1, ub1, ...		% lower and upper bounds for the 1st layer
	lb2, ub2, ...		% lower and upper bounds for the 2nd layer
	lb3, ub3, ...		% lower and upper bounds for the 3rd layer
	options1, ...		% options for the 1st layer
	options2, ...		% options for the 2nd layer
	options3)			% options for the 3rd layer
% MMMDEB1B_C Sequential Max-Min-Max DE/best/1/bin with finding constraint
% boundary 
% [xbest1, xbest2, xbest3] = MMMDEB1B_C(fitfun, maxfunevals, lb1,
% ub1, lb2, ub2, lb3, ub3) maximizes the function fitfun associated with a
% maximizer xbest1 among box limitations [lb1, ub1], a minimizer xbest2
% among [lb2, ub2], and a maximizer xbest3 among [lb3, ub3], which are
% searched by evolutionary algorithm within a maximal function evaluations
% maxfunevals.
% MMMDEB1B_C(..., options1) maximizes the function with the given
% options Options1 for the 1st layer.
% MMMDEB1B_C(..., options1, options2) maximizes the function with
% the given options Options2 for the 2nd layer.
% MMMDEB1B_C(..., options1, options2, options3) maximizes the
% function with the given options Options3 for the 3rd layer.
% [..., fbest] = MMMDEB1B_C(...) returns the function value of the
% max-min-max solution.
if nargin <= 9
	options1 = [];
end

if nargin <= 10
	options2 = [];
end

if nargin <= 11
	options3 = [];
end

D1 = numel(lb1);
D2 = numel(lb2);
D3 = numel(lb3);

% Default options for Layer 1
defaultOptions1.dimensionFactor = 10;
defaultOptions1.F = 0.9;
defaultOptions1.CR = 0.5;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolX = 0;
defaultOptions1.TolFun = 0;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.InnerSolver = 'mmdeb1b_pce';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.innerState = [];

defaultOptions1.TolCon = 1e-6;
defaultOptions1.nonlcon = [];
defaultOptions1.innerMaxIter = 200;
defaultOptions1.migrateFactor = 0.6;
defaultOptions1.zeta = 1e-4;

options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2.dimensionFactor = 10;
defaultOptions2.F = 0.9;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.TolFun = 0;
defaultOptions2.TolX = 0;
defaultOptions2.InnerSolver = 'debest1bin';
options2 = setdefoptions(options2, defaultOptions2);

% Default options for Layer 3
defaultOptions3.dimensionFactor = 10;
options3 = setdefoptions(options3, defaultOptions3);

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
zeta = options1.zeta;

X = options1.initial.X;
f = options1.initial.f;
innerState = options1.initial.innerState;
existInnerState = ~isempty(innerState);

NP1 = ceil(dimensionFactor * D1);
NP2 = ceil(options2.dimensionFactor * D2);
NP3 = ceil(options3.dimensionFactor * D3);

% Initialize contour data
if isDisplayIter
	plotFitfun = @(x, y) feval(fitfun, x, y, 0.5 * (lb3 + ub3));
	[XX, YY, ZZ] = minmaxcontourdata(D1, lb1, ub1, lb2, ub2, plotFitfun);
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
innerXbest1 = zeros(D2, NP1);
innerXbest2 = zeros(D3, NP1);
innerUbest1 = innerXbest1;
innerUbest2 = innerXbest2;
innerXCbest1 = innerXbest1;
innerUCbest1 = innerXbest1;
innerXCbest2 = innerXbest2;
innerUCbest2 = innerXbest2;
V = X;
U = X;
PX2 = zeros(D2, NP2, NP1);
PX3 = zeros(D3, NP3, NP2, NP1);
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
	innerMaxfunevalsX = innerMaxIter * NP2 * NP3;
	
	for i = 1 : NP1
		innerFitfun = @(y, z) feval(fitfun, X(:, i), y, z);	
		optionsX2i = options2;
		
		if ~isempty(nonlcon)
			innerNonlcon = @(y, z) feval(nonlcon, X(:, i), y, z);
			optionsX2i.nonlcon = innerNonlcon;
		end
		
		if existInnerState
			optionsX2i.initial = innerState{i};			
		end
		
		[innerXbest1(:, i), innerXbest2(:, i), innerFbest, innerOutX{i}] = ...
			feval(innerSolver, innerFitfun, innerMaxfunevalsX, lb2, ub2, ...
			lb3, ub3, optionsX2i, options3);
		
		f(i) = -innerFbest;
		innerState{i} = innerOutX{i}.final;
	end
	
	for i = 1 : NP1		
		counteval = counteval + innerOutX{i}.fes(end);
		countcon = countcon + innerOutX{i}.countcon;
	end
end

% Constraint violation measure
cm = zeros(1, NP1);
nc = zeros(1, NP1);

for i = 1 : NP1
	clb = lb1 - X(:, i);
	cub = X(:, i) - ub1;
	cm(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
	nc(i) = sum(clb > 0) + sum(cub > 0);
end

for i = 1 : NP1
	clb = lb2 - innerXbest1(:, i);
	cub = innerXbest1(:, i) - ub2;
	cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
	nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
end

for i = 1 : NP1
	clb = lb3 - innerXbest2(:, i);
	cub = innerXbest2(:, i) - ub3;
	cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
	nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
end

if ~isempty(nonlcon)
	for i = 1 : NP1
		[cx, ceqx] = ...
			feval(nonlcon, X(:, i), innerXbest1(:, i), innerXbest2(:, i));
		
		countcon = countcon + 1;		
		cm(i) = cm(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
		nc(i) = nc(i) + sum(cx > 0) + sum(ceqx > 0);
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
innerXbest1 = innerXbest1(:, pfidx);
innerXbest2 = innerXbest2(:, pfidx);
innerState = innerState(pfidx);
cm = cm(pfidx);
nc = nc(pfidx);

% Display
if isDisplayIter
	if all(isinf(f))
		displayitermessages([X; innerXbest1], [U; innerUbest1], ...
			cm, countiter, ...
			XX, YY, ZZ, 'counteval', counteval, ...
			'successRate', successRate);
	else
		displayitermessages([X; innerXbest1], [U; innerUbest1], ...
			f(~isinf(f)), countiter, ...
			XX, YY, ZZ, 'counteval', counteval, ...
			'successRate', successRate);
	end
	
	display_inner_info(innerState);
end

% Record minimal function values
out = updateoutput(out, X, f, counteval, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs));

countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval >= maxfunevals;
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
		for j = 1 : NP2
			if isempty(innerState{i}.innerState{j})
				anyEmptyInnerState = true;
				break;
			end
		end
		
		if anyEmptyInnerState || isempty(innerState{i})
			anyEmptyInnerState = true;
			break;
		end
	end
	
	if ~anyEmptyInnerState		
		for i = 1 : NP1			
			% Copy from itselfs individuals
			for j = 1 : NP2
				PX2(:, j, i) = innerState{i}.X(:, j);
				
				for k = 1 : NP3
					PX3(:, k, j, i) = innerState{i}.innerState{j}.X(:, k);
				end
			end
			
			% Copy from innerXbest
			migrationNP2 = ceil(migrateFactor * NP2);
			beginIndex1 = NP2 - migrationNP2;
			for j = 1 : migrationNP2
				r = floor(NP1 * rand + 1);
				PX2(:, beginIndex1 + j, i) = innerXbest1(:, r) .* (1 + zeta * randn(D2, 1));
			end
			
			migrationNP3 = ceil(migrateFactor * NP3);
			beginIndex2 = NP3 - migrationNP3;
			for j = 1 : NP2
				for k = 1 : migrationNP3
					r = floor(NP1 * rand + 1);
					PX3(:, beginIndex2 + k, j, i) = innerXbest2(:, r) .* (1 + zeta * randn(D3, 1));
				end
			end
		end
	else		
		for i = 1 : NP1
			for j = 1 : NP2
				PX2(:, j, i) = lb2 + (ub2 - lb2) .* rand(D2, 1);
				
				for k = 1 : NP3
					PX3(:, k, j, i) = lb3 + (ub3 - lb3) .* rand(D3, 1);
				end
			end
		end
	end
	
	PU2 = PX2;
	PU3 = PX3;
	
	if ~isempty(nonlcon)
		parfor i = 1 : NP1
			% Compute XC
			innerFitfunXCi = @(y, z) max(feval(nonlcon, X(:, i), y , z)).^2;
			
			innerOptionsXCi = options2;
			innerOptionsXCi.initial = [];
			
			[innerXCbest1(:, i), innerXCbest2(:, i), ~, innerOutXC{i}] = ...
				feval(innerSolver, innerFitfunXCi, innerMaxfunevalsX, ...
				lb2, ub2, ...
				lb3, ub3, ...
				innerOptionsXCi, options3);
			
			% Compute UC
			innerFitfunUCi = @(y, z) max(feval(nonlcon, U(:, i), y, z)).^2;
			
			innerOptionsUCi = options2;
			innerOptionsUCi.initial = [];
			
			[innerUCbest1(:, i), innerUCbest2(:, i), ~, innerOutUC{i}] = ...
				feval(innerSolver, innerFitfunUCi, innerMaxfunevalsX, ...
				lb2, ub2, ...
				lb3, ub3, ...
				innerOptionsXCi, options3);
		end
		
		for i = 1 : NP1
			PX2(:, end, i) = innerXCbest1(:, i);
			PU2(:, end, i) = innerUCbest1(:, i);
			
			for j = 1 : NP2
				PX3(:, end, j, i) = innerXCbest2(:, i);
				PU3(:, end, j, i) = innerUCbest2(:, i);
			end
			
			countcon = countcon + innerOutXC{i}.fes(end);
			countcon = countcon + innerOutUC{i}.fes(end);
		end
	end
	
	% Selection
	innerMaxfunevalsX = innerMaxIter * NP2 * NP3;
	for i = 1 : NP1
		% Compute fxi, f(i)
		innerFitfunXi = @(y, z) feval(fitfun, X(:, i), y, z);
		optionsX2i = options2;
		optionsX2i.initial = [];		
		optionsX2i.initial.X = PX2(:, :, i);
		
		for j = 1 : NP2
			optionsX2i.initial.innerState{j} = [];
			optionsX2i.initial.innerState{j}.X = PX3(:, :, j, i);
		end
		
		if ~isempty(nonlcon)
			optionsX2i.nonlcon = @(y, z) feval(nonlcon, X(:, i), y, z);
		end
		
		[innerXbest1(:, i), innerXbest2(:, i), innerFbest, innerOutX{i}] = ...
			feval(innerSolver, innerFitfunXi, innerMaxfunevalsX, ...
			lb2, ub2, ...
			lb3, ub3, ...
			optionsX2i, options3);
		
		X_Converged_FEs(i) = innerOutX{i}.fes(end);
		f(i) = -innerFbest;
		
		% Compute fui
		innerFitfunUi = @(y, z) feval(fitfun, U(:, i), y, z);
		optionsU2i = options2;
		optionsU2i.initial = [];
		optionsU2i.initial.X = PU2(:, :, i);
		
		for j = 1 : NP2
			optionsU2i.initial.innerState{j} = [];
			optionsU2i.initial.innerState{j}.X = PU3(:, :, j, i);
		end
		
		if ~isempty(nonlcon)
			optionsU2i.nonlcon = @(y, z) feval(nonlcon, U(:, i), y, z);
		end
		
		[innerUbest1(:, i), innerUbest2(:, i), innerFbest, innerOutU{i}] = ...
			feval(innerSolver, innerFitfunUi, innerMaxfunevalsX, ...
			lb2, ub2, ...
			lb3, ub3, ...
			optionsU2i, options3);
		
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
		clb = lb2 - innerXbest1(:, i);
		cub = innerXbest1(:, i) - ub2;
		cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
				
		clb = lb2 - innerUbest1(:, i);
		cub = innerUbest1(:, i) - ub2;
		cm_u(i) = cm_u(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = nc_u(i) + sum(clb > 0) + sum(cub > 0);
	end
	
	for i = 1 : NP1
		clb = lb3 - innerXbest2(:, i);
		cub = innerXbest2(:, i) - ub3;
		cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
				
		clb = lb3 - innerUbest2(:, i);
		cub = innerUbest2(:, i) - ub3;
		cm_u(i) = cm_u(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = nc_u(i) + sum(clb > 0) + sum(cub > 0);
	end
	
	if ~isempty(nonlcon)
		for i = 1 : NP1
			[cx, ceqx] = feval(...
				nonlcon, ...
				X(:, i), ...
				innerXbest1(:, i), ...
				innerXbest2(:, i));
			
			countcon = countcon + 1;
			cm(i) = cm(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
			nc(i) = nc(i) + sum(cx > 0) + sum(ceqx > 0);
			
			[cu, cequ] = feval(...
				nonlcon, ...
				U(:, i), ...
				innerUbest1(:, i), ...
				innerUbest2(:, i));
			
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
			innerXbest1(:, i) = innerUbest1(:, i);
			innerXbest2(:, i) = innerUbest2(:, i);
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
			displayitermessages([X; innerXbest1], [U; innerUbest1], ...
				cm, countiter, ...
				XX, YY, ZZ, 'counteval', counteval, ...
				'successRate', successRate);
		else
			displayitermessages([X; innerXbest1], [U; innerUbest1], ...
				f(~isinf(f)), countiter, ...
				XX, YY, ZZ, 'counteval', counteval, ...
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
	innerXbest1 = innerXbest1(:, pfidx);
	innerXbest2 = innerXbest2(:, pfidx);
	innerState = innerState(pfidx);
	
	% Record
	out = updateoutput(out, X, f, counteval, ...
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
xbest3 = innerState{1}.innerState{1}.X(:, 1);
fbest = -f(1);

final.innerState = innerState;

out = finishoutput(out, X, f, counteval, 'final', final, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'countcon', countcon);
end