function [xmin, fmin, out] = dcmaeabin(fitfun, lb, ub, maxfunevals, options)
% DCMAEABIN Differential covariance matrix adaptation evolutionary
% algorithm
% DCMAEABIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DCMAEABIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 1;
defaultOptions.mu_CR = 0.9;
defaultOptions.Display = 'off';
defaultOptions.FactorNP = 1;
defaultOptions.Restart = 0;
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.X = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
mu_CR = defaultOptions.mu_CR;
isDisplayIter = strcmp(options.Display, 'iter');
FactorNP = max(1, options.FactorNP);
Restart = max(0, round(options.Restart));
RecordPoint = max(1, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
X = options.X;

D = numel(lb);
chiN = D^0.5 * (1 - 1 / (4 * D) + 1 / (21 * D^2));

if isempty(X)
	NP = ceil(dimensionFactor * D);
else
	[~, NP] = size(X);
end

% Too small maxfunevals Exception
if maxfunevals < NP
	NP = maxfunevals;
	X = zeros(D, NP);
	f = zeros(1, NP);
	if maxfunevals < 1e3
		LHS = lhsdesign(maxfunevals, D, 'iteration', 50)';
	elseif maxfunevals < 1e4
		LHS = lhsdesign(maxfunevals, D, 'iteration', 5)';
	else
		LHS = rand(D, maxfunevals);
	end
	
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* LHS(:, i);
		f(i) = feval(fitfun, X(:, i));
	end
	
	% Prepare output arguments
	[fmin, fminidx] = min(f);
	xmin = X(:, fminidx);
	
	% Display
	if isDisplayIter
		countiter = 1;
		displayitermessages(X, f, countiter);
	end
		
	% Output info
	out = initoutput(RecordPoint, D, NP, maxfunevals);
	out = finishoutput(out, X, f, maxfunevals);
	out.bestever.fmin = fmin;
	out.bestever.xmin = xmin;
	return;
end

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D, NP, maxfunevals);
initNP = NP;
MAX_NP = 2000;

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = preparecontourdata(D, lb, ub, fitfun);
end

for iRestart = 1 : (Restart + 1)
	NP = min(MAX_NP, round(initNP * FactorNP ^ (iRestart - 1)));
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
% 	sigma = 0.5;
% 	diagD = (ub - lb) / max(ub - lb);
	diagD = ones(D, 1);
	diagC = diagD.^2;
	C = diag(diagC);
	B = eye(D, D);
	BD = B .* repmat(diagD', D, 1);
	m = (lb + ub) / 2;
% 	m = rand(D, 1);
	Z = randn(D, NP);

	% Initialize population
	if iRestart >= 2 || isempty(X)
		X = zeros(D, NP);
		for i = 1 : NP
			X(:, i) = m + sigma * BD * Z(:, i);
		end
	end
	
	% Repair
	for i = 1 : NP
		for j = 1 : D			
			while X(j, i) < lb(j) || X(j, i) > ub(j)
				if X(j, i) < lb(j)
					X(j, i) = 2 * lb(j) - X(j, i);
				elseif X(j, i) > ub(j)
					X(j, i) = 2 * ub(j) - X(j, i);
				end
			end
		end
	end
	
	% Evaluation
	f = zeros(1, NP);
	counteval = counteval + NP;
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
	end
	
	% Initialize variables
	V = X;
	U = X;
	
	% Display
	if isDisplayIter
		displayitermessages(X, U, f, countiter, XX, YY, ZZ);
	end
	
	% Record
	out = updateoutput(out, X, f, counteval);
	
	% Iteration counter
	countiter = countiter + 1;
	
	while true		
		% Termination conditions
		stdf = std(f);
		stdX = std(X, 0, 2);
		meanX = mean(X, 2);
		outofmaxfunevals = counteval > maxfunevals - NP;
		reachftarget = min(f) <= ftarget;
		fitnessconvergence = stdf <= mean(f) * 100 * eps || stdf <= realmin || stdf <= TolFun;
		solutionconvergence = all(stdX <= meanX * 100 * eps) || all(stdX <= 100 * realmin) || ...
			all(stdX <= TolX);
		
		% Convergence conditions
		if outofmaxfunevals || reachftarget || fitnessconvergence || solutionconvergence
			break;
		end
		
		% Sort population
		[f, findex] = sort(f);
		X = X(:, findex);
		Z = Z(:, findex);
	
		% Update the variables of CMA-ES	
		zmean = Z(:, 1 : mu) * w;
		mold = m;
		m = X(:, 1 : mu) * w;
		ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * (B * zmean);
% 		hsig = norm(ps) / sqrt(1 - (1 - cs)^(2 * countiter)) / chiN < 1.4 + 2 / (D + 1);
% 		pc = (1 - cc) * pc ...
% 			+ hsig * (sqrt(cc * (2 - cc) * mueff) / sigma) * (m - mold);
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
% 		P = 0.2;
		F = 0.5 + 0.5 * rand(1, NP);
		
		% Mutation
		for i = 1 : NP
			r1 = floor(1 + NP * rand);
			
			while r1 == i
				r1 = floor(1 + NP * rand);
			end
			
			r2 = floor(1 + NP * rand);
			
			while r2 == i || r2 == r1
				r2 = floor(1 + NP * rand);
			end
			
			Z(:, i) = randn(D, 1);
			V(:, i) = m * (1 - P) + P * X(:, i) + P * F(i) * (X(:, r1) - X(:, r2)) + ...
				(1 - P) * sigma * BD * Z(:, i);
		end
		
		for i = 1 : NP
			% Binominal Crossover
			jrand = floor(1 + D * rand);
			CR = mu_CR + 0.1 * randn;
			
			while CR > 1
				CR = mu_CR + 0.1 * randn;
			end
				
			for j = 1 : D
				
				if rand < CR || j == jrand
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
					U(j, i) = X(j, i) + rand * (lb(j) - X(j, i));
				elseif U(j, i) > ub(j)
					U(j, i) = X(j, i) + rand * (ub(j) - X(j, i));
				end
			end
		end
		
		% Display
		if isDisplayIter
			displayitermessages(X, U, f, countiter, XX, YY, ZZ);
		end
		
		% Selection		
		counteval = counteval + NP;
		for i = 1 : NP
			fui = feval(fitfun, U(:, i));
			
			if fui < f(i)
				f(i) = fui;
				X(:, i) = U(:, i);
			end
		end
		
		% Record
		out = updateoutput(out, X, f, counteval);
		
		% Iteration counter
		countiter = countiter + 1;
	end
	
	[fmin, fminidx] = min(f);
	xmin = X(:, fminidx);
	
	if fmin < out.bestever.fmin
		out.bestever.fmin = fmin;
		out.bestever.xmin = xmin;
	end
	
	% Termination conditions
	if counteval > maxfunevals - NP
		break;
	end
	
	if min(f) <= ftarget
		break;
	end
end

out = finishoutput(out, X, f, counteval);
fmin = out.bestever.fmin;
xmin = out.bestever.xmin;
end