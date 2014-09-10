function [xbest1, xbest2, fbest, out] = mmdeb1b_c(fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, options1, options2)
% MMDEB1B_C Sequential Min-Max DE/best/1/bin with constraint activation
% MMDEB1B_C(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MMDEB1B_C(..., options1) minimizes the function with the given
% options Options1 for the 1st layer.
% MMDEB1B_C(..., options1, options2) minimize the function with the
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
defaultOptions1.NP = 10;
defaultOptions1.F = 0.9;
defaultOptions1.CR = 0.5;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.InnerSolver = 'debest1bin';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.psai = [];
defaultOptions1.initial.innerState = [];

defaultOptions1.nonlcon = [];
defaultOptions1.innerMaxIter = 200;
defaultOptions1.migrateFactor = 0.7;
defaultOptions1.ConstraintHandling = 'EpsilonMethod';
defaultOptions1.EpsilonValue = 0;
defaultOptions1.EarlyStop = 'auto';

options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2.NP = 10;
defaultOptions2.F = 0.9;
defaultOptions2.CR = 0.5;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.ConstraintHandling = 'EpsilonMethod';
defaultOptions2.EpsilonValue = 0;
defaultOptions2.EarlyStop = 'auto';
options2 = setdefoptions(options2, defaultOptions2);

% Initialize algorithmic variables
F = options1.F;
CR = options1.CR;
isDisplayIter = strcmp(options1.Display, 'iter');
RecordPoint = max(0, floor(options1.RecordPoint));
TolStagnationIteration = options1.TolStagnationIteration;
innerSolver = options1.InnerSolver;
nonlcon = options1.nonlcon;
innerMaxIter = options1.innerMaxIter;
migrateFactor = options1.migrateFactor;

EpsilonValue = options1.EpsilonValue;
if ~isempty(strfind(options1.ConstraintHandling, 'EpsilonMethod'))
	EpsilonMethod = true;
else
	EpsilonMethod = false;
end

if ~isempty(strfind(options1.EarlyStop, 'auto'))
	EarlyStop = true;
else
	EarlyStop = false;
end

if ~isempty(options1.initial)
	options1.initial = setdefoptions(options1.initial, defaultOptions1.initial);
	X = options1.initial.X;
	fx = options1.initial.f;
	innerState = options1.initial.innerState;
	psai_x = options1.initial.psai;
else
	X = [];
	fx = [];
	innerState = [];
	psai_x = [];
end

existInnerState = ~isempty(innerState);

NP1 = options1.NP;
NP2 = options2.NP;

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
psai_u = zeros(1, NP1);

out = initoutput(RecordPoint, D1, NP1, maxfunevals, ...
	'innerFstd', ...
	'innerMeanXstd', ...
	'successRate', ...
	'X_Converged_FEs', ...
	'U_Converged_FEs');

% Evaluation
if isempty(fx)
	fx = zeros(1, NP1);
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
		
		fx(i) = -innerFbest;
		innerState{i} = innerOutX{i}.final;
	end
	
	for i = 1 : NP1		
		counteval = counteval + innerOutX{i}.fes(end);
		countcon = countcon + innerOutX{i}.countcon;
	end
end

% Constraint violation
if isempty(psai_x) && EpsilonMethod
	psai_x = zeros(1, NP1);
	for i = 1 : NP1		
		clbx1 = lb1 - X(:, i);
		cubx1 = X(:, i) - ub1;
		psai_x(i) = sum(clbx1(clbx1 > 0)) + sum(cubx1(cubx1 > 0));
				
		clbx2 = lb2 - innerXbest(:, i);
		cubx2 = innerXbest(:, i) - ub2;
		psai_x(i) = psai_x(i) + sum(clbx2(clbx2 > 0)) + sum(cubx2(cubx2 > 0));
		
		if ~isempty(nonlcon)			
			[cx, ceqx] = feval(nonlcon, X(:, i), innerXbest(:, i));
			countcon = countcon + 1;
			psai_x(i) = psai_x(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
		end
	end
end

% Sort
if ~EpsilonMethod
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	innerXbest = innerXbest(:, fidx);
	innerState = innerState(fidx);
else
	PsaiFx = [psai_x', fx'];
	[~, SortingIndex] = sortrows(PsaiFx);
	X = X(:, SortingIndex);
	fx = fx(SortingIndex);
	innerXbest = innerXbest(:, SortingIndex);
	innerState = innerState(SortingIndex);
	psai_x = psai_x(SortingIndex);
end

% Display
if isDisplayIter
	if all(isinf(fx))
		dispconitermsg([X; innerXbest], [U; innerUbest], ...
			psai_x, countiter, ...
			XX, YY, ZZ, CC, 'counteval', counteval, ...
			'successRate', successRate);
	else
		dispconitermsg([X; innerXbest], [U; innerUbest], ...
			fx(~isinf(fx)), countiter, ...
			XX, YY, ZZ, CC, 'counteval', counteval, ...
			'successRate', successRate);
	end
	
	display_inner_info(innerState);
end

% Record minimal function values
out = updateoutput(out, X, fx, counteval + countcon, countiter, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs));

countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval + countcon >= maxfunevals;
	if ~EarlyStop
		if outofmaxfunevals
			break;
		end
	else		
		TolX = 10 * eps(mean(X(:)));
		solutionconvergence = std(X(:)) <= TolX;
		TolFun = 10 * eps(mean(fx));
		functionvalueconvergence = std(fx(:)) <= TolFun;
		stagnation = countStagnation >= TolStagnationIteration;
		
		if outofmaxfunevals || ...
				solutionconvergence || ...
				functionvalueconvergence || ...
				stagnation
			break;
		end
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
% 			n_migration = ceil(migrateFactor * NP2);
% 			beginIndex = NP2 - n_migration;
% 			for j = 1 : n_migration
% 				r = floor(NP1 * rand + 1);
% 				PX2(:, beginIndex + j, i) = innerXbest(:, r);
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
% 		fx(i) = -innerFbest;
		
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
% 		counteval = counteval + innerOutX{i}.fes(end);
% 		countcon = countcon + innerOutX{i}.countcon;
		counteval = counteval + innerOutU{i}.fes(end);
		countcon = countcon + innerOutU{i}.countcon;
	end

	% Constraint violation
	if EpsilonMethod
		for i = 1 : NP1
			clbu1 = lb1 - U(:, i);
			cubu1 = U(:, i) - ub1;
			psai_u(i) = sum(clbu1(clbu1 > 0)) + sum(cubu1(cubu1 > 0));
			
			clbu2 = lb2 - innerUbest(:, i);
			cubu2 = innerUbest(:, i) - ub2;
			psai_u(i) = psai_u(i) + sum(clbu2(clbu2 > 0)) + sum(cubu2(cubu2 > 0));
			
			if ~isempty(nonlcon)
				[cu, cequ] = feval(nonlcon, U(:, i), innerUbest(:, i));
				countcon = countcon + 1;
				psai_u(i) = psai_u(i) + sum(cu(cu > 0)) + sum(cequ(cequ > 0));
			end
		end
	end
	
	% Replacement
	successRate = 0;
	FailedIteration = true;
	if ~EpsilonMethod
		for i = 1 : NP1			
			if fu(i) < fx(i)
				fx(i) = fu(i);
				X(:, i) = U(:, i);
				innerXbest(:, i) = innerUbest(:, i);
				innerState{i} = innerOutU{i}.final;
				successRate = successRate + 1 / NP1;
				FailedIteration = false;
			else
				innerState{i} = innerOutX{i}.final;
			end
		end
	else
		for i = 1 : NP1
			X_AND_U_IN_EPSILON = psai_u(i) < EpsilonValue && psai_x(i) < EpsilonValue;
			X_AND_U_EQUAL_EPSILON = psai_u(i) == psai_x(i);
			
			if ((X_AND_U_IN_EPSILON || X_AND_U_EQUAL_EPSILON) && fu(i) < fx(i)) || ...
					(~X_AND_U_IN_EPSILON && psai_u(i) < psai_x(i))
				fx(i) = fu(i);
				X(:, i) = U(:, i);
				innerXbest(:, i) = innerUbest(:, i);
				innerState{i} = innerOutU{i}.final;
				successRate = successRate + 1 / NP1;
				FailedIteration = false;
				psai_x(i)	= psai_u(i);
			else
				innerState{i} = innerOutX{i}.final;
			end
		end
	end
	
	% Display
	if isDisplayIter
		if all(isinf(fx))
			dispconitermsg([X; innerXbest], [U; innerUbest], ...
				psai_x, countiter, ...
				XX, YY, ZZ, CC, 'counteval', counteval, ...
				'successRate', successRate);
		else
			dispconitermsg([X; innerXbest], [U; innerUbest], ...
				fx(~isinf(fx)), countiter, ...
				XX, YY, ZZ, CC, 'counteval', counteval, ...
				'successRate', successRate);
		end
		
		display_inner_info(innerState);
	end
	
	% Sort
	if ~EpsilonMethod
		[fx, fidx] = sort(fx);
		X = X(:, fidx);
		innerXbest = innerXbest(:, fidx);
		innerState = innerState(fidx);
	else
		PsaiFx = [psai_x', fx'];
		[~, SortingIndex] = sortrows(PsaiFx);
		X = X(:, SortingIndex);
		fx = fx(SortingIndex);
		innerXbest = innerXbest(:, SortingIndex);
		innerState = innerState(SortingIndex);
		psai_x = psai_x(SortingIndex);
	end
	
	% Record
	out = updateoutput(out, X, fx, counteval + countcon, countiter, ...
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
fbest = fx(1);

final.innerState = innerState;

out = finishoutput(out, X, fx, counteval + countcon, countiter, ...
	'final', final, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'countcon', countcon);
end