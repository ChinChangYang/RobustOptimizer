function [xmin, fmin, out] = rbdebin(fitfun, lb, ub, maxfunevals, options)
% RBDEBIN Differential evolution with rank-based mutation
% RBDEBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% RBDEBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.beta = 3.0;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.ConstraintHandling = 'Interpolation';

options = setdefoptions(options, defaultOptions);
F = options.F;
CR = options.CR;
beta = options.beta;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

if isequal(options.ConstraintHandling, 'Interpolation')
	interpolation = true;
else
	interpolation = false;
end

D = numel(lb);

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	fx = options.initial.f;
else
	X = [];
	fx = [];
end

if isempty(X)	
	NP = options.NP;
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'FC');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
end

% Initialize population
if isempty(X)
	X = zeros(D, NP);
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* rand(D, 1);
	end
end

% Evaluation
if isempty(fx)
	fx = zeros(1, NP);
	for i = 1 : NP
		fx(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
end

% Sort
[fx, fidx] = sort(fx);
X = X(:, fidx);

% Initialize variables
V = X;
U = X;
fu = zeros(1, NP);
FC = zeros(1, NP);		% Consecutive Failure Counter
rt = zeros(1, NP);
r1 = zeros(1, NP);
r2 = zeros(1, NP);
r3 = zeros(1, NP);

% Display
if isDisplayIter
	displayitermessages(...
		X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter, ...
	'FC', FC);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(fx) <= ftarget;
	stagnation = countStagnation >= TolStagnationIteration;	
	if outofmaxfunevals || reachftarget || stagnation
		break;
	end
	
	% SV-based framework
	for i = 1 : NP
		rt(i) = i;
		
		% Generate r1
		r1(i) = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		while rt(i) == r1(i)
			r1(i) = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		end
		
		% Generate r2
		r2(i) = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		while rt(i) == r2(i) || r1(i) == r2(i)
			r2(i) = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		end
		
		% Generate r3
		r3(i) = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		while rt(i) == r3(i) || r1(i) == r3(i) || r2(i) == r3(i)
			r3(i) = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		end
	end
	
	% Mutation
	for i = 1 : NP		
		V(:, i) = X(:, r1(i)) + F * (X(:, r2(i)) - X(:, r3(i)));
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		for j = 1 : D
			if rand < CR || j == jrand
				U(j, i) = V(j, i);
			else
				U(j, i) = X(j, rt(i));
			end
		end
	end
		
	if interpolation
		% Correction for outside of boundaries
		for i = 1 : NP
			for j = 1 : D
				if U(j, i) < lb(j)
					U(j, i) = 0.5 * (lb(j) + X(j, rt(i)));
				elseif U(j, i) > ub(j)
					U(j, i) = 0.5 * (ub(j) + X(j, rt(i)));
				end
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			X, U, fx, countiter, XX, YY, ZZ);
	end
	
	% Evaluation
	for i = 1 : NP
		fu(i) = feval(fitfun, U(:, i));
		counteval = counteval + 1;
	end
	
	% Selection
	FailedIteration = true;
	for i = 1 : NP		
		if fu(i) < fx(i)
			X(:, i)		= U(:, i);
			fx(i)		= fu(i);
			FailedIteration = false;
			FC(i)		= 0;
		else
			FC(i) = FC(i) + 1;
		end
	end
	
	% Sort
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	FC = FC(fidx);
	
	% Record
	out = updateoutput(out, X, fx, counteval, countiter, ...
		'FC', FC);
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end	
end

[fmin, minindex] = min(fx);
xmin = X(:, minindex);

out = finishoutput(out, X, fx, counteval, countiter, ...
	'FC', zeros(NP, 1));
end