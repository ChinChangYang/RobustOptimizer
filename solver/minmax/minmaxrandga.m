function [xminmax1, xminmax2, fminmax, out] = minmaxrandga(fitfun1, ...
	maxfunevals1, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXRANDGA Classical DE/rand/1/bin for min-max problems
% MINMAXRANDGA(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% MINMAXRANDGA(fitfun, lb, ub, maxfunevals, D_Min, D_Max) requires that
% the dimension D_Min and D_Max of variables for minimizer and maximizer,
% respectively.
% MINMAXRANDGA(..., options) minimize the function by solver options.
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
D = D_X + D_Y;
NP_X = max(ceil(dimensionFactor * D_X), floor(maxfunevalsFactor * maxfunevals));
NP_Y = 2 * max(ceil(dimensionFactor * D_Y), floor(maxfunevalsFactor * maxfunevals));
NP_S = NP_X * NP_Y;
lb_X = lb(1 : D_X);
lb_Y = lb(D_X + 1 : D_X + D_Y);
ub_X = ub(1 : D_X);
ub_Y = ub(D_X + 1 : D_X + D_Y);
lb_S = lb;
ub_S = ub;

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
YU = Y;
realSamples = 1;
fw = zeros(2 * NP_X, 2 * NP_Y);

% Evaluation
f = zeros(NP_X, NP_Y);
counteval = counteval + NP_S;
for i = 1 : NP_X
	for j = 1 : NP_Y
		f(i, j) = feval(fitfun, X(:, i), Y(:, j));
	end
end

[fmax, fmaxidx] = max(f, [], 2);
[fminmax, fminidx] = min(fmax);

% Display
if isDisplayIter
	S = zeros(D_X + D_Y, NP_S);
	SU = zeros(D_X + D_Y, NP_S);
	for i = 1 : NP_X
		for j = 1 : NP_Y
			S(:, (i - 1)* NP_Y + j) = [X(:, i); Y(:, j)];
			SU(:, (i - 1)* NP_Y + j) = [XU(:, i); YU(:, j)];
		end
	end
	displayitermessages(S, SU, f(:), countiter, XX, YY, ZZ);
end

% Record minimal function values
out = updateoutput(out, X, fmax, counteval);
countiter = countiter + 1;

while true	
	% Termination conditions
	stdf = std(f(:));
	stdX1 = std(X, 0, 2);
	meanX1 = mean(X, 2);
	stdX2 = std(Y, 0, 2);
	meanX2 = mean(Y, 2);
	outofmaxfunevals = counteval > maxfunevals - NP_S;
	fitnessconvergence = stdf <= mean(f(:)) * 100 * eps || stdf <= realmin || stdf <= TolFun;
	solutionconvergence1 = all(stdX1 <= meanX1 * 100 * eps) || all(stdX1 <= 100 * realmin) || ...
		all(stdX1 <= TolX);
	solutionconvergence2 = all(stdX2 <= meanX2 * 100 * eps) || all(stdX2 <= 100 * realmin) || ...
		all(stdX2 <= TolY);
	
	if outofmaxfunevals || fitnessconvergence || ...
			solutionconvergence1 || solutionconvergence2
		break;
	end
		
	% X's Mutation
	for i = 1 : NP_X
		r1 = floor(1 + NP_X * rand);
		r2 = floor(1 + NP_X * rand);
		r3 = r1;
		
		while r1 == r3 || r2 == r3
			r3 = floor(1 + NP_X * rand);
		end
		
		XV(:, i) = X(:, r1) + (F + 0.05 * randn) * (X(:, r2) - X(:, r3));
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
	
	% Randomize Y
	for i = 1 : NP_Y
		YU(:, i) = lb_Y + (ub_Y - lb_Y) .* rand(D_Y, 1);
	end
	YU(:, 1) = lb_Y;
	YU(:, 2) = ub_Y;
	
	% Display
	if isDisplayIter
		S = zeros(D_X + D_Y, NP_S);
		SU = zeros(D_X + D_Y, NP_S);
		for i = 1 : NP_X
			for j = 1 : NP_Y
				S(:, (i - 1)* NP_Y + j) = [X(:, i); Y(:, j)];
				SU(:, (i - 1)* NP_Y + j) = [XU(:, i); YU(:, j)];
			end
		end
		displayitermessages(S, SU, f(:), countiter, XX, YY, ZZ);
	end
	
	% Selection
	XW = [X, XU];
	YW = [Y, YU];
	counteval = counteval + 4 * NP_S;
	for i = 1 : 2 * NP_X
		for j = 1 : 2 * NP_Y
			fw(i, j) = feval(fitfun, XW(:, i), YW(:, j));
		end
	end
		
	[fwmax, fmaxidx] = max(fw, [], 2);
	
	for i = 1 : NP_X
		if fwmax(NP_X+i) < fwmax(i)
			X(:, i) = XW(:, NP_X+i);
		end
	end
	
	Y = YU;
	
	[fminmax, fminidx] = min(fwmax);
	
	% Record
	out = updateoutput(out, XW, fwmax, counteval);
	countiter = countiter + 1;
end

out = finishoutput(out, XW, fwmax, counteval);
xminmax1 = XW(:, fminidx);
xminmax2 = YW(:, fmaxidx(fminidx));
end