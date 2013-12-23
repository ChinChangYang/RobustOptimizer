function [xmin, fmin, out] = degleig(fitfun, lb, ub, maxfunevals, options)
% DEGLEIG DEGL with eigenvector-based crossover
% DEGLEIG(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DEGLEIG(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 10;
defaultOptions.maxfunevalsFactor =  0;
defaultOptions.R = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.F = 0.7;
defaultOptions.NeighborhoodRatio = 0.1;
defaultOptions.Display = 'off';
defaultOptions.FactorNP = 1;
defaultOptions.Restart = 0;
defaultOptions.RecordPoint = 100;
defaultOptions.Noise = false;
defaultOptions.SampleFactor = 1.01;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = 0;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.X = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
R = options.R;
CR = options.CR;
F = options.F;
isDisplayIter = strcmp(options.Display, 'iter');
FactorNP = max(1, options.FactorNP);
Restart = max(0, round(options.Restart));
RecordPoint = max(1, floor(options.RecordPoint));
Noise = options.Noise;
SampleFactor = options.SampleFactor;
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
alpha = F;
beta = F;
X = options.X;

D = numel(lb);

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
	
	% Initialize population
	if iRestart >= 2 || isempty(X)
		if NP < 1e1
			LHS = lhsdesign(NP, D, 'iteration', 10)';
		elseif NP < 1e2
			LHS = lhsdesign(NP, D, 'iteration', 2)';
		else
			LHS = rand(D, NP);
		end
		
		X = zeros(D, NP);
		for i = 1 : NP
			X(:, i) = lb + (ub - lb) .* LHS(:, i);
		end
	end
	
	% Initialize variables
	k = ceil(0.5 * (options.NeighborhoodRatio * NP));
	w = 0.05 + 0.9 * rand(1, NP);
	wc = w;
	V = X;
	U = X;
	XT = X;
	VT = X;
	UT = X;
	realSamples = 1;
	
	% Evaluation
	f = zeros(1, NP);
	counteval = counteval + NP;
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
	end
	
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
		fitnessconvergence = stdf <= mean(abs(f)) * 100 * eps || stdf <= realmin || stdf <= TolFun;
		solutionconvergence = all(stdX <= meanX * 100 * eps) || all(stdX <= 100 * realmin) || ...
			all(stdX <= TolX);
		
		% Convergence conditions
		if outofmaxfunevals || reachftarget || fitnessconvergence || solutionconvergence
			break;
		end
		
		% Mutation
		% Global best
		[~, g_best] = min(f);
		
		for i = 1 : NP
			% Neiborhoods index
			n_index = (i-k) : (i+k);
			lessthanone = n_index < 1;
			n_index(lessthanone) = n_index(lessthanone) + NP;
			greaterthanNP = n_index > NP;
			n_index(greaterthanNP) = n_index(greaterthanNP) - NP;
			
			% Neiborhood solutions and fitness
			Xn = X(:, n_index);
			fn = f(n_index);
			
			% Best neiborhood
			[~, n_besti] = min(fn);
			Xn_besti = Xn(:, n_besti);
			
			% Random neiborhood index
			n_index(n_index == i) = [];
			Xn = X(:, n_index);
			p = ceil(rand * numel(n_index));			
			q = ceil(rand * numel(n_index));
			
			while p == q				
				q = ceil(rand * numel(n_index));
			end
			
			% Random neiborhood solutions
			Xp = Xn(:, p);
			Xq = Xn(:, q);
			
			% Local donor vector
			Li = X(:, i) + alpha * (Xn_besti - X(:, i)) + ...
				beta * (Xp - Xq);
			
			% Global donor vector
			r1 = floor(1 + NP * rand);
			
			while i == r1
				r1 = floor(1 + NP * rand);
			end
			
			r2 = floor(1 + NP * rand);
			
			while i == r2 || r1 == r2
				r2 = floor(1 + NP * rand);
			end
			
			gi = X(:, i) + alpha * (X(:, g_best) - X(:, i)) + ...
				beta * (X(:, r1) - X(:, r2));
						
			% Self-adaptive weight factor
			wc(i) = w(i) + F * (w(g_best) - w(i)) + ...
				F * (w(r1) - w(r2));
			
			if wc(i) < 0.05
				wc(i) = 0.05;
			elseif wc(i) > 0.95
				wc(i) = 0.95;
			end
			
			V(:, i) = wc(i) * gi + (1 - wc(i)) * Li;
		end
		
		[B, ~] = eig(cov(X'));
		for i = 1 : NP
			if rand < R
				% Rotational Crossover
				XT(:, i) = B' * X(:, i);
				VT(:, i) = B' * V(:, i);
				jrand = floor(1 + D * rand);
				
				for j = 1 : D
					if rand < CR || j == jrand
						UT(j, i) = VT(j, i);
					else
						UT(j, i) = XT(j, i);
					end
				end
				
				U(:, i) = B * UT(:, i);
			else
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
		if Noise
			prevSamples = ceil(realSamples);
			realSamples = SampleFactor * realSamples;
			samples = ceil(realSamples);
			for i = 1 : NP
				if samples - prevSamples >= 1
					% Compute fxi, f(i)
					fxi = 0.0;
					for j = 1 : samples - prevSamples
						fxi = fxi + feval(fitfun, X(:, i));
						counteval = counteval + 1;
					end
					fxi = fxi / (samples - prevSamples);
					f(i) = (prevSamples / samples) * f(i) + ...
						(1 - prevSamples / samples) * fxi;
				end
				
				% Compute fui
				fui = 0.0;
				for j = 1 : samples
					fui = fui + feval(fitfun, U(:, i));
					counteval = counteval + 1;
				end
				fui = fui / samples;
				
				% Replacement
				if fui < f(i)
					f(i) = fui;
					X(:, i) = U(:, i);
					w(i) = wc(i);
				end
			end
		else
			counteval = counteval + NP;
			for i = 1 : NP
				fui = feval(fitfun, U(:, i));
				
				if fui < f(i)
					f(i) = fui;
					X(:, i) = U(:, i);
					w(i) = wc(i);
				end
			end			
		end
		
		% Record
		out = updateoutput(out, X, f, counteval);
		
		% Iteration counter
		countiter = countiter + 1;
	end
	
	[fmin, fminidx] = min(f);
	xmin = X(:, fminidx);
	
	if Noise
		out.bestever.fmin = mean(f);
		out.bestever.xmin = mean(X, 2);
	else
		if fmin < out.bestever.fmin
			out.bestever.fmin = fmin;
			out.bestever.xmin = xmin;
		end
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