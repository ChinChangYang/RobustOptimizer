function [xmin, fmin, out] = rbde_sps(fitfun, lb, ub, maxfunevals, options)
% RBDE_SPS Differential evolution with rank-based mutation and SPS
% Framework 
% RBDE_SPS(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% RBDE_SPS(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.Q = 70;
defaultOptions.beta = 3.0;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];

options = setdefoptions(options, defaultOptions);
F = options.F;
CR = options.CR;
beta = options.beta;
Q = options.Q;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

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
SP = X;
fSP = fx;
iSP = 1;
[~, sortidxfSP] = sort(fSP);

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
		
	% Mutation
	for i = 1 : NP
		% Generate r1
		r1 = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		while i == r1
			r1 = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		end
		
		% Generate r2
		r2 = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		while i == r2 || r1 == r2
			r2 = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		end
		
		% Generate r3
		r3 = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		while i == r3 || r1 == r3 || r2 == r3
			r3 = floor(1 + NP/2/(beta-1)*(beta-sqrt(beta^2-4*(beta-1)*rand)));
		end
		
		if FC(i) <= Q
			V(:, i) = X(:, r1) + F * (X(:, r2) - X(:, r3));
		else
			V(:, i) = SP(:, sortidxfSP(r1)) + ...
				F * (SP(:, sortidxfSP(r2)) - SP(:, sortidxfSP(r3)));			
		end
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		
		if FC(i) <= Q
			for j = 1 : D
				if rand < CR || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, i);
				end
			end
		else
			for j = 1 : D
				if rand < CR || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = SP(j, i);
				end
			end
		end
	end
	
	% Correction for outside of boundaries
	for i = 1 : NP
		if FC(i) <= Q
			for j = 1 : D
				if U(j, i) < lb(j)
					U(j, i) = 0.5 * (lb(j) + X(j, i));
				elseif U(j, i) > ub(j)
					U(j, i) = 0.5 * (ub(j) + X(j, i));
				end
			end
		else
			for j = 1 : D
				if U(j, i) < lb(j)
					U(j, i) = 0.5 * (lb(j) + SP(j, i));
				elseif U(j, i) > ub(j)
					U(j, i) = 0.5 * (ub(j) + SP(j, i));
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
			SP(:, iSP)	= U(:, i);
			fSP(iSP)	= fu(i);
			iSP			= mod(iSP, NP) + 1;
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
	[~, sortidxfSP] = sort(fSP);
	
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