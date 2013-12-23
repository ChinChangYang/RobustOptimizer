function [xminmax1, xminmax2, fminmax, out] = minmaxjade(fitfun1, ...
	maxfunevals1, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXJADE Two-level JADE for min-max problems
% MINMAXJADE(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% MINMAXJADE(fitfun, lb, ub, maxfunevals, D_Min, D_Max) requires that
% the dimension D_Min and D_Max of variables for minimizer and maximizer,
% respectively.
% MINMAXJADE(..., options) minimize the function by solver options.
if nargin <= 6
	options = [];
else
	options = options1;
end

fitfun = fitfun1;
maxfunevals = maxfunevals1;
lb = [lb1; lb2];
ub = [ub1; ub2];
D_Min = numel(lb1);
D_Max = numel(lb2);

defaultOptions.dimensionFactor = 5;
defaultOptions.maxfunevalsFactor = 0;
defaultOptions.mu_CR = 0.5;
defaultOptions.mu_F = 1.0;
defaultOptions.delta_CR = 0.1;
defaultOptions.delta_F = 0.1;
defaultOptions.p = 0.05;
defaultOptions.w = 0.1;
defaultOptions.Display = 'off';
defaultOptions.FactorNP = 1;
defaultOptions.Restart = 0;
defaultOptions.RecordPoint = 100;
defaultOptions.Noise = false;
defaultOptions.SampleFactor = 1.01;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;

options = setdefoptions(options, defaultOptions);

options_Y = options;
options_Y.Display = 'off';
options_Y.TolX = 1e-8;
options_Y.TolFun = 1e-8;
maxfunevals_Y = maxfunevals;

dimensionFactor = max(1, options.dimensionFactor);
maxfunevalsFactor = options.maxfunevalsFactor;
mu_CR = options.mu_CR;
mu_F = options.mu_F;
delta_CR = options.delta_CR;
delta_F = options.delta_F;
p = options.p;
w = options.w;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
TolFun = options.TolFun;
TolX = options.TolX;

D = D_Min;
NP = max(ceil(dimensionFactor * D), floor(maxfunevalsFactor * maxfunevals));
lb_Y = lb(D_Min + 1 : D_Min + D_Max);
ub_Y = ub(D_Min + 1 : D_Min + D_Max);
lb = lb(1 : D_Min);
ub = ub(1 : D_Min);

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D_Min, NP, maxfunevals1);

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = minmaxcontourdata(numel(lb1), lb1, ub1, lb2, ub2, fitfun1);
end

% Initialize population
X = zeros(D, NP);
for i = 1 : NP
	X(:, i) = lb + (ub - lb) .* rand(D, 1);
end

% Initialize variables
V = X;
U = X;
Y = zeros(D_Max, NP);
YU = Y;
pbest_size = p * NP;

% Evaluation
f = zeros(1, NP);
for i = 1 : NP
	nest_fitfun = @(y) -feval(fitfun, X(:, i), y);
	[Y(:, i), fmax, nestout] = jadebin(nest_fitfun, lb_Y, ub_Y, maxfunevals_Y, options_Y);
	counteval = counteval + nestout.fes(end);
	f(i) = -fmax;
end

[f, fidx] = sort(f);
X = X(:, fidx);

% Initialize archive
A = zeros(D, 2 * NP);
A(:, 1 : NP) = X;

% Display
if isDisplayIter
	displayitermessages([X; Y], [U; Y], f, countiter, XX, YY, ZZ, 'mu_F', mu_F, 'mu_CR', mu_CR);
end

% Record minimal function values
out = updateoutput(out, X, fmax, counteval);
countiter = countiter + 1;

while true
	% Termination conditions
	stdf = std(f);
	stdX1 = std(X, 0, 2);
	meanX1 = mean(X, 2);
	outofmaxfunevals = counteval > maxfunevals1 - NP;
	fitnessconvergence = stdf <= mean(f) * 100 * eps || stdf <= realmin || stdf <= TolFun;
	solutionconvergence = all(stdX1 <= meanX1 * 100 * eps) || all(stdX1 <= 100 * realmin) || ...
		all(stdX1 <= TolX);
	
	% Convergence conditions
	if outofmaxfunevals || fitnessconvergence || solutionconvergence
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
	
	F(F <= 0) = eps;
	
	A_Counter = 0;
	XA = [X, A];
	
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
			
			V(:, i) = X(:, i) + F(i) * (X(:, pbest_idx) - X(:, i) + X(:, r1) - XA(:, r2));
			
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
	
	% Repair
	for i = 1 : NP
		for j = 1 : D
			if U(j, i) < lb(j)
				U(j, i) = 2 * lb(j) - U(j, i);
			elseif U(j, i) > ub(j)
				U(j, i) = 2 * ub(j) - U(j, i);
			else
				continue;
			end
			
			if U(j, i) < lb(j)
				U(j, i) = lb(j);
			elseif U(j, i) > ub(j)
				U(j, i) = ub(j);
			end
		end
	end
	
	% Selection
	for i = 1 : NP
		nest_fitfun = @(y) -feval(fitfun, U(:, i), y);
		[YU(:, i), fmax, nestout] = jadebin(nest_fitfun, lb_Y, ub_Y, maxfunevals_Y, options_Y);
		counteval = counteval + nestout.fes(end);
		fui = -fmax;
		
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
			Y(:, i) = YU(:, i);
			A(:, NP + A_Counter + 1) = U(:, i);
			S_CR(A_Counter + 1) = CR(i);
			S_F(A_Counter + 1) = F(i);
			A_Counter = A_Counter + 1;
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			[X; Y], [U; YU], f, countiter, XX, YY, ZZ, ...
			'mu_F', mu_F, 'mu_CR', mu_CR);
	end
	
	% Update archive
	rand_idx = randperm(NP + A_Counter);
	A(:, 1 : NP) = A(:, rand_idx(1 : NP));
	
	% Update CR and F
	if A_Counter > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : A_Counter));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : A_Counter).^2) / sum(S_F(1 : A_Counter));
	end
	
	% Sort
	[f, fidx] = sort(f);
	X = X(:, fidx);
	
	% Record
	out = updateoutput(out, X, f, counteval);	
	countiter = countiter + 1;
end

out = finishoutput(out, X, f, counteval);
xminmax1 = X(:, 1);
xminmax2 = Y(:, 1);
fminmax = min(f);
end