function [xminmax1, xminmax2, fminmax, out] = minmaxdebest1bin(fitfun1, ...
	maxfunevals1, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXDEBEST1BIN Bi-level DE/best/1/bin for min-max problems
% MINMAXDEBEST1BIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% MINMAXDEBEST1BIN(fitfun, lb, ub, maxfunevals, D_Min, D_Max) requires that
% the dimension D_Min and D_Max of variables for minimizer and maximizer,
% respectively.
% MINMAXDEBEST1BIN(..., options) minimize the function by solver options.
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
defaultOptions.maxfunevalsFactor =  0;
defaultOptions.CR = 0.9;
defaultOptions.F = 0.7;
defaultOptions.Display = 'off';
defaultOptions.FactorNP = 1;
defaultOptions.Restart = 0;
defaultOptions.RecordPoint = 100;
defaultOptions.Noise = true;
defaultOptions.SampleFactor = 1.04;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;

options = setdefoptions(options, defaultOptions);
options_Y = options;
options_Y.dimensionFactor = 5;
options_Y.CR = 0.9;
options_Y.F = 0.7;
options_Y.Display = 'off';
options_Y.Noise = false;
options_Y.TolX = 1e-8;
options_Y.TolFun = 1e-8;
maxfunevals_Y = 100 * D_Max;

dimensionFactor = options.dimensionFactor;
maxfunevalsFactor = options.maxfunevalsFactor;
CR = options.CR;
F = options.F;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
Noise = options.Noise;
SampleFactor = options.SampleFactor;
TolFun = options.TolFun;
TolX = options.TolX;

D = D_Min;
NP = max(ceil(dimensionFactor * D), floor(maxfunevalsFactor * maxfunevals));
lb_Y = lb(D_Min + 1 : D_Min + D_Max);
ub_Y = ub(D_Min + 1 : D_Min + D_Max);
lb_S = lb;
ub_S = ub;
lb = lb(1 : D_Min);
ub = ub(1 : D_Min);

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D, NP, maxfunevals);

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = minmaxcontourdata(numel(lb1), lb1, ub1, lb2, ub2, fitfun1);
end

X = zeros(D, NP);
for i = 1 : NP
	X(:, i) = lb + (ub - lb) .* rand(D, 1);
end

% Initialize variables
V = X;
U = X;
Y = zeros(D_Max, NP);
YU = Y;

% Evaluation
f = zeros(1, NP);
for i = 1 : NP
	nest_fitfun = @(y) -feval(fitfun, X(:, i), y);
	[Y(:, i), fmax, nestout] = debest1bin(nest_fitfun, lb_Y, ub_Y, maxfunevals_Y, options_Y);
	counteval = counteval + nestout.fes(end);
	f(i) = -fmax;
end

% Display
if isDisplayIter
	displayitermessages([X; Y], [U; Y], f, countiter, XX, YY, ZZ);
end

out = updateoutput(out, X, f, counteval);
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
	
	% Mutation
	[~, best] = min(f);
	for i = 1 : NP
		r1 = floor(1 + NP * rand);
		r2 = r1;
		
		while r1 == r2
			r2 = floor(1 + NP * rand);
		end
		
		V(:, i) = X(:, best) + (F + 0.05 * randn) * (X(:, r1) - X(:, r2));
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
	
	% Repair X and U
	for i = 1 : NP
		for j = 1 : D
			while X(j, i) < lb(j)
				X(j, i) = 2 * lb(j) - X(j, i);
			end
			
			while X(j, i) > ub(j)
				X(j, i) = 2 * ub(j) - X(j, i);
			end
			
			while U(j, i) < lb(j)
				U(j, i) = 2 * lb(j) - U(j, i);
			end
			
			while U(j, i) > ub(j)
				U(j, i) = 2 * ub(j) - U(j, i);
			end
		end
	end
	
	% Prepare YR
	% 		YR = zeros(D_Max, dimensionFactor * D_Max);
	% 		for i = 1 : dimensionFactor * D_Max
	% 			YR(:, i) = lb_Y + (ub_Y - lb_Y) * rand(D_Max, 1);
	% 		end
	% 		options_Y.InitialPopulation = [Y(:, best), YR];
	
	% Selection
	if Noise
		maxfunevals_Y = floor(SampleFactor * maxfunevals_Y) + D_Max;
		for i = 1 : NP
			% Compute fxi, f(i)
			nest_fitfun = @(y) -feval(fitfun, X(:, i), y);
			[Y(:, i), fmax, nestout] = debest1bin(nest_fitfun, lb_Y, ub_Y, maxfunevals_Y, options_Y);
			counteval = counteval + nestout.fes(end);
			f(i) = -fmax;
			
			% Compute fui
			nest_fitfun = @(y) -feval(fitfun, U(:, i), y);
			[YU(:, i), fmax, nestout] = debest1bin(nest_fitfun, lb_Y, ub_Y, maxfunevals_Y, options_Y);
			counteval = counteval + nestout.fes(end);
			fui = -fmax;
			
			% Replacement
			if fui < f(i)
				f(i) = fui;
				X(:, i) = U(:, i);
				Y(:, i) = YU(:, i);
			end
		end
	else
		for i = 1 : NP
			nest_fitfun = @(y) -feval(fitfun, U(:, i), y);
			[YU(:, i), fmax, nestout] = debest1bin(nest_fitfun, lb_Y, ub_Y, maxfunevals, options_Y);
			counteval = counteval + nestout.fes(end);
			fui = -fmax;
			
			if fui < f(i)
				f(i) = fui;
				X(:, i) = U(:, i);
				Y(:, i) = YU(:, i);
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages([X; Y], [U; YU], f, countiter, XX, YY, ZZ);
	end
	
	out = updateoutput(out, X, f, counteval);
	countiter = countiter + 1;
end

[fminmax, fminidx] = min(f);
xminmax1 = X(:, fminidx);
xminmax2 = Y(:, fminidx);

out = finishoutput(out, X, f, counteval);
end