function [xminmax1, xminmax2, fminmax, out] = minmaxcodebest1bin(fitfun1, ...
	maxfunevals1, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXDERAND1BIN Co-evol. DE/best/1/bin for min-max problems
% MINMAXDERAND1BIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% MINMAXDERAND1BIN(fitfun, lb, ub, maxfunevals, D_Min, D_Max) requires that
% the dimension D_Min and D_Max of variables for minimizer and maximizer,
% respectively.
% MINMAXDERAND1BIN(..., options) minimize the function by solver options.
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
defaultOptions.CR = 0.5;
defaultOptions.F = 0.9;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
maxfunevalsFactor = options.maxfunevalsFactor;
CR = options.CR;
F = options.F;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
TolFun = options.TolFun;
TolX = options.TolX;
TolY = TolX;

D_X = D_Min;
D_Y = D_Max;
NP_X = max(ceil(dimensionFactor * D_X), floor(maxfunevalsFactor * maxfunevals));
NP_Y = max(ceil(dimensionFactor * D_Y), floor(maxfunevalsFactor * maxfunevals));
NP_S = NP_X * NP_Y;
lb_X = lb(1 : D_X);
lb_Y = lb(D_X + 1 : D_X + D_Y);
ub_X = ub(1 : D_X);
ub_Y = ub(D_X + 1 : D_X + D_Y);

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D_Min, NP_X, maxfunevals1);

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = minmaxcontourdata(numel(lb1), lb1, ub1, lb2, ub2, fitfun1);
end

% Initialize population
if NP_X < 1e1
	LHS_X = lhsdesign(NP_X, D_X, 'iteration', 10)';
elseif NP_X < 1e2
	LHS_X = lhsdesign(NP_X, D_X, 'iteration', 2)';
else
	LHS_X = rand(D_X, NP_X);
end

if NP_Y < 1e1
	LHS_Y = lhsdesign(NP_Y, D_Y, 'iteration', 10)';
elseif NP_Y < 1e2
	LHS_Y = lhsdesign(NP_Y, D_Y, 'iteration', 2)';
else
	LHS_Y = rand(D_Y, NP_Y);
end

X = zeros(D_X, NP_X);
Y = zeros(D_Y, NP_Y);
for i = 1 : NP_X
	X(:, i) = lb_X + (ub_X - lb_X) .* LHS_X(:, i);
end
for i = 1 : NP_Y
	Y(:, i) = lb_Y + (ub_Y - lb_Y) .* LHS_Y(:, i);
end

% Initialize variables
XV = X;
XU = X;
YV = Y;
YU = Y;
fu = zeros(NP_X, NP_Y);

% Evaluation
f = zeros(NP_X, NP_Y);
counteval = counteval + NP_S;
for i = 1 : NP_X
	for j = 1 : NP_Y
		f(i, j) = feval(fitfun, X(:, i), Y(:, j));
	end
end

[fmin, ~] = min(f);
[fmax, fmaxidx] = max(f, [], 2);
[fminmax, fminidx] = min(fmax);

% Display
if isDisplayIter
	S = zeros(D_X + D_Y, NP_X * NP_Y);
	for i = 1 : NP_X
		for j = 1 : NP_Y
			S(:, (i - 1)* NP_Y + j) = [X(:, i); Y(:, j)];
		end
	end
	displayitermessages(S, S, f, countiter, XX, YY, ZZ);
end

out = updateoutput(out, X, fmax, counteval);
countiter = countiter + 1;

while true
	% Termination conditions
	stdfmax = std(fmax);
	stdfmin = std(fmin);
	stdX1 = std(X, 0, 2);
	meanX1 = mean(X, 2);
	stdX2 = std(Y, 0, 2);
	meanX2 = mean(Y, 2);
	outofmaxfunevals = counteval > maxfunevals - NP_X * NP_Y;
	fitnessconvergence1 = stdfmax <= mean(fmax) * 100 * eps || stdfmax <= realmin || stdfmax <= TolFun;
	fitnessconvergence2 = stdfmin <= mean(fmin) * 100 * eps || stdfmin <= realmin || stdfmin <= TolFun;
	solutionconvergence1 = all(stdX1 <= meanX1 * 100 * eps) || all(stdX1 <= 100 * realmin) || ...
		all(stdX1 <= TolX);
	solutionconvergence2 = all(stdX2 <= meanX2 * 100 * eps) || all(stdX2 <= 100 * realmin) || ...
		all(stdX2 <= TolY);
	
	if outofmaxfunevals || fitnessconvergence1 || fitnessconvergence2 || ...
			solutionconvergence1 || solutionconvergence2
		break;
	end
	
	% X's Mutation
	[~, fminmaxindex] = min(fmax);
	for i = 1 : NP_X
		r1 = floor(1 + NP_X * rand);
		r2 = floor(1 + NP_X * rand);
		
		while r1 == r2
			r2 = floor(1 + NP_X * rand);
		end
		
		XV(:, i) = X(:, fminmaxindex) + (F + 0.01 * randn) * (X(:, r1) - X(:, r2));
	end
	
	% X's Crossover
	for i = 1 : NP_X
		jrand = floor(1 + D_X * rand);
		for j = 1 : D_X
			if rand < CR || j == jrand
				XU(j, i) = XV(j, i);
			else
				XU(j, i) = X(j, i);
			end
		end
	end
	
	% Y's Mutation
	[~, fmaxminindex] = max(fmin);
	for i = 1 : NP_Y
		r1 = floor(1 + NP_Y * rand);
		r2 = floor(1 + NP_Y * rand);
		
		while r1 == r2
			r2 = floor(1 + NP_Y * rand);
		end
		
		YV(:, i) = Y(:, fmaxminindex) + (F + 0.01 * randn) * (Y(:, r1) - Y(:, r2));
	end
	
	% Y's Crossover
	for i = 1 : NP_Y		
		jrand = floor(1 + D_Y * rand);
		for j = 1 : D_Y
			if rand < CR || j == jrand
				YU(j, i) = YV(j, i);
			else
				YU(j, i) = Y(j, i);
			end
		end
	end
	
	% X's Repair
	for i = 1 : NP_X
		for j = 1 : D_X
			if XU(j, i) < lb_X(j)
				XU(j, i) = 2 * lb_X(j) - XU(j, i);
			elseif XU(j, i) > ub_X(j)
				XU(j, i) = 2 * ub_X(j) - XU(j, i);
			else
				continue;
			end
			
			if XU(j, i) < lb_X(j)
				XU(j, i) = lb_X(j);
			elseif XU(j, i) > ub_X(j)
				XU(j, i) = ub_X(j);
			end
		end
	end
	
	% Y's Repair
	for i = 1 : NP_Y
		for j = 1 : D_Y
			if YU(j, i) < lb_Y(j)
				YU(j, i) = 2 * lb_Y(j) - YU(j, i);
			elseif YU(j, i) > ub_Y(j)
				YU(j, i) = 2 * ub_Y(j) - YU(j, i);
			else
				continue;
			end
			
			if YU(j, i) < lb_Y(j)
				YU(j, i) = lb_Y(j);
			elseif YU(j, i) > ub_Y(j)
				YU(j, i) = ub_Y(j);
			end
		end
	end
	
	% Display
	if isDisplayIter
		S = zeros(D_X + D_Y, NP_X * NP_Y);
		for i = 1 : NP_X
			for j = 1 : NP_Y
				S(:, (i - 1)* NP_Y + j) = [X(:, i); Y(:, j)];
			end
		end
		displayitermessages(S, S, f, countiter, XX, YY, ZZ);
	end
	
	% Selection
	counteval = counteval + NP_S;
	for i = 1 : NP_X
		for j = 1 : NP_Y
			fu(i, j) = feval(fitfun, XU(:, i), YU(:, j));
		end
	end
		
	[fumax, ~] = max(fu, [], 2);
	
	for i = 1 : NP_X
		if fumax(i) < fmax(i)
			X(:, i) = XU(:, i);
		end
	end
	
	[fumin, ~] = min(fu);
	
	for i = 1 : NP_Y
		if fumin(i) > fmin(i)
			Y(:, i) = YU(:, i);
		end
	end
	
	counteval = counteval + NP_S;
	for i = 1 : NP_X
		for j = 1 : NP_Y
			f(i, j) = feval(fitfun, X(:, i), Y(:, j));
		end
	end
	
	[fmin, ~] = min(f);
	[fmax, fmaxidx] = max(f, [], 2);
	[fminmax, fminidx] = min(fmax);
	
	out = updateoutput(out, X, fmax, counteval);
	
	countiter = countiter + 1;
end

xminmax1 = X(:, fminidx);
xminmax2 = Y(:, fmaxidx(fminidx));
out = finishoutput(out, X, fmax, counteval);
end