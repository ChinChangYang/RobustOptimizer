function [xmin, fmin, out] = dcmaea_sps(fitfun, lb, ub, maxfunevals, options)
% DCMAEA_SPS DCMA-EA with SPS Framework
% DCMAEA_SPS(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DCMAEA_SPS(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.CR = 0.5;
defaultOptions.Q = 70;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.ConstraintHandling = 'Interpolation';

options = setdefoptions(options, defaultOptions);
mu_CR = options.CR;
Q = options.Q;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolStagnationIteration = options.TolStagnationIteration;

if isequal(options.ConstraintHandling, 'Interpolation')
	interpolation = true;
else
	interpolation = false;
end

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	X = options.initial.X;
	fx = options.initial.f;
else
	X = [];
	fx = [];
end

D = numel(lb);
chiN = D^0.5 * (1 - 1 / (4 * D) + 1 / (21 * D^2));
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
mu = floor(0.5 * NP);
w = log(mu + 0.5)-log(1:mu)';
mueff = sum(w)^2 / sum(w.^2);
w = w / sum(w);
cs = (mueff + 2) / (D + mueff + 5);
ds = 1 + 2 * max(0, sqrt((mueff - 1) / (D + 1)) - 1) + cs;
cc = (4 + mueff / D) / (D + 4 + 2 * mueff / D);
c1 = 2 / ((D + 1.3)^2 + mueff);
cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((D + 2)^2 + mueff));
ps = zeros(D, 1);
pc = zeros(D, 1);
sigma = max(ub - lb) / 4;
diagD = ones(D, 1);
diagC = diagD.^2;
C = diag(diagC);
B = eye(D, D);
% BD = B .* repmat(diagD', D, 1);
m = (lb + ub) / 2;
Z = randn(D, NP);

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
	
	% Update the variables of CMA-ES
	zmean = Z(:, 1 : mu) * w;
	mold = m;
	m = X(:, 1 : mu) * w;
	ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * (B * zmean);
	pc = (1 - cc) * pc + sqrt(cc * (2 - cc) * mueff) * (m - mold) / sigma;
	arpos = (X(:, 1 : mu) - repmat(mold, 1, mu)) / sigma;
	C = (1 - c1 - cmu) * C ...
		+ c1 * (pc * pc') ...
		+ cmu ...
		* arpos * (repmat(w, 1, D) .* arpos');
	sigma = sigma * exp((sqrt(sum(ps.^2)) / chiN - 1) * cs / ds);
	C = triu(C) + triu(C,1)';
	C = real(C);
	[B, temp] = eig(C);
	B = real(B);
	temp = real(temp);
	diagD = diag(temp);
	
	% Covariance matrix repair
	if min(diagD) <= 0
		diagD(diagD < 0) = 0;
		temp = 1e-14 * max(diagD);
		C = C + temp * eye(D, D);
		diagD = diagD + temp * ones(D,1);
	end
	
	if max(diagD) > 1e14 * min(diagD)
		temp = 1e-14 * max(diagD) - min(diagD);
		C = C + temp * eye(D, D);
		diagD = diagD + temp * ones(D,1);
	end
	
	diagD = sqrt(diagD);
	BD = B .* repmat(diagD', D, 1);
	
	% Update the variables of DE
	P = 0.5 * (1 + counteval / maxfunevals);
	F = 0.5 + 0.5 * rand(1, NP);
		
	% Crossover rates
	CR = mu_CR + 0.1 * randn(1, NP);
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	
	% Mutation
	for i = 1 : NP		
		% Generate r1
		r1 = floor(1 + NP * rand);
		while i == r1
			r1 = floor(1 + NP * rand);
		end
		
		% Generate r2
		r2 = floor(1 + NP * rand);
		while i == r2 || r1 == r2
			r2 = floor(1 + NP * rand);
		end
		
		Z(:, i) = randn(D, 1);
		
		if FC(i) <= Q
			V(:, i) = m * (1 - P) + P * X(:, i) + ...
				P * F(i) * (X(:, r1) - X(:, r2)) + ...
				(1 - P) * sigma * BD * Z(:, i);
		else
			V(:, i) = m * (1 - P) + P * SP(:, i) + ...
				P * F(i) * (SP(:, r1) - SP(:, r2)) + ...
				(1 - P) * sigma * BD * Z(:, i);
		end
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		if FC(i) <= Q
			for j = 1 : D
				if rand < CR(i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = X(j, i);
				end
			end
		else
			for j = 1 : D
				if rand < CR(i) || j == jrand
					U(j, i) = V(j, i);
				else
					U(j, i) = SP(j, i);
				end
			end			
		end
	end
	
	if interpolation
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
			FC(i)		= 0;
			FailedIteration = false;
		else
			FC(i) = FC(i) + 1;
		end
	end
	
	% Sort	
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	Z = Z(:, fidx);
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

[fmin, fminidx] = min(fx);
xmin = X(:, fminidx);

out = finishoutput(out, X, fx, counteval, countiter, ...
	'FC', zeros(NP, 1));
end