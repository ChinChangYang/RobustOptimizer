function [xmin, fmin, out] = sadeeig(fitfun, lb, ub, maxfunevals, options)
% SADEEIG SaDE/eig
% SADEEIG(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% SADEEIG(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.maxfunevalsFactor = 0;
defaultOptions.R = 0.5;
defaultOptions.LP = 40;
defaultOptions.Display = 'off';
defaultOptions.FactorNP = 1;
defaultOptions.Restart = 0;
defaultOptions.RecordPoint = 100;
defaultOptions.Noise = false;
defaultOptions.SampleFactor = 1.01;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.X = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
R = options.R;
LP = round(options.LP);
isDisplayIter = strcmp(options.Display, 'iter');
FactorNP = max(1, options.FactorNP);
Restart = max(0, round(options.Restart));
RecordPoint = max(1, floor(options.RecordPoint));
Noise = options.Noise;
SampleFactor = options.SampleFactor;
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
D = numel(lb);
X = options.X;

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
	V = X;
	U = X;
	XT = X;
	VT = X;
	UT = X;
	realSamples = 1;	
	K = 4;
	ns = ones(K, LP);
	nf = ones(K, LP);
	CRMemory = cell(1, K);
	for i = 1 : K
		CRMemory{i} = 0.5 + 0.1 * randn;
	end
	selStrategy = zeros(1, NP);
	CR = zeros(K, NP);
	CRm = 0.5 * ones(K, 1);
	p = 1/K * ones(K, 1);
	g = 2;
	
	% Evaluation
	f = zeros(1, NP);
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
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
		
		if countiter > LP
			S = sum(ns, 2) ./ (sum(ns, 2) + sum(nf, 2)) + eps;
			p = S ./ sum(S);
			ns(:, 1:end-1) = ns(:, 2:end);
			nf(:, 1:end-1) = nf(:, 2:end);
		end
		
		cp = zeros(K, 1);
		cp(1) = p(1);
		for i = 2 : K
			cp(i) = cp(i-1) + p(i);
		end
		
		for i = 1 : NP
			s = rand;
			selStrategy(i) = find(cp >= s, 1, 'first');
		end
		
		F = 0.5 + 0.3 * randn(1, NP);
		
		if g >= LP
			for k = 1 : K
				CRm(k) = median(CRMemory{k});
			end
		end
		
		for k = 1 : K
			CR(k, :) = CRm(k) + 0.1 * randn(1, NP);
			checkIndex = CR(k, :) < 0 | CR(k, :) > 1;
			while any(checkIndex > 0)
				CR(k, checkIndex) = CRm(k) + 0.1 * randn(1, sum(checkIndex));
				checkIndex = CR(k, :) < 0 | CR(k, :) > 1;
			end
		end
		
		% Mutation
		[~, besti] = min(f);
		for i = 1 : NP
			if selStrategy(i) == 1
				% DE/rand/1
				r1 = floor(1 + NP * rand);
				r2 = floor(1 + NP * rand);
				r3 = r2;
				
				while r2 == r3
					r3 = floor(1 + NP * rand);
				end
				
				pF = 1 + 0.01 * randn;
				V(:, i) = X(:, r1) + F(i) * pF * (X(:, r2) - X(:, r3));
			elseif selStrategy(i) == 2
				% DE/current-to-best/2
				r2 = floor(1 + NP * rand);
				r3 = r2;
				r4 = floor(1 + NP * rand);
				r5 = r4;
				
				while r2 == r3
					r3 = floor(1 + NP * rand);
				end
				
				while r4 == r5
					r5 = floor(1 + NP * rand);
				end
				
				pF1 = 1 + 0.03 * randn;
				pF2 = 1 + 0.03 * randn;
				pF3 = 1 + 0.03 * randn;
				V(:, i) = X(:, i) + F(i) * pF1 * (X(:, besti) - X(:, i)) + ...
					F(i) * pF2 * (X(:, r2) - X(:, r3)) + ...
					F(i) * pF3 * (X(:, r4) - X(:, r5));
			elseif selStrategy(i) == 3
				% DE/rand/2
				r1 = floor(1 + NP * rand);
				r2 = floor(1 + NP * rand);
				r3 = r2;
				r4 = floor(1 + NP * rand);
				r5 = r4;
				
				while r2 == r3
					r3 = floor(1 + NP * rand);
				end
				
				while r4 == r5
					r5 = floor(1 + NP * rand);
				end
				
				pF1 = 1 + 0.03 * randn;
				pF2 = 1 + 0.03 * randn;
				V(:, i) = X(:, r1) + F(i) * pF1 * (X(:, r2) - X(:, r3)) + ...
					F(i) * pF2 * (X(:, r4) - X(:, r5));
			elseif selStrategy(i) == 4
				% DE/current-to-rand/1
				r1 = floor(1 + NP * rand);
				r2 = floor(1 + NP * rand);
				r3 = r2;
				
				while r2 == r3
					r3 = floor(1 + NP * rand);
				end
				
				pF1 = 1 + 0.03 * randn;
				pF2 = 1 + 0.03 * randn;
				V(:, i) = X(:, i) + rand * pF1 * (X(:, r1) - X(:, i)) + ...
					F(i) * pF2 * (X(:, r2) - X(:, r3));
			else
				fprintf('BUG\n');
			end
		end
		
		[B, ~] = eig(cov(X'));
		for i = 1 : NP
			if selStrategy(i) == 1
				% DE/rand/1/bin
				if rand < R
					% Rotational Crossover
					XT(:, i) = B' * X(:, i);
					VT(:, i) = B' * V(:, i);
					jrand = floor(1 + D * rand);
					
					for j = 1 : D
						if rand < CR(1, i) || j == jrand
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
						if rand < CR(1, i) || j == jrand
							U(j, i) = V(j, i);
						else
							U(j, i) = X(j, i);
						end
					end
				end
			elseif selStrategy(i) == 2
				% DE/rand-to-best/1/bin
				if rand < R
					% Rotational Crossover
					XT(:, i) = B' * X(:, i);
					VT(:, i) = B' * V(:, i);
					jrand = floor(1 + D * rand);
					
					for j = 1 : D
						if rand < CR(2, i) || j == jrand
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
						if rand < CR(2, i) || j == jrand
							U(j, i) = V(j, i);
						else
							U(j, i) = X(j, i);
						end
					end
				end
			elseif selStrategy(i) == 3
				% DE/rand/2/bin
				if rand < R
					% Rotational Crossover
					XT(:, i) = B' * X(:, i);
					VT(:, i) = B' * V(:, i);
					jrand = floor(1 + D * rand);
					
					for j = 1 : D
						if rand < CR(3, i) || j == jrand
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
						if rand < CR(3, i) || j == jrand
							U(j, i) = V(j, i);
						else
							U(j, i) = X(j, i);
						end
					end
				end
			elseif selStrategy(i) == 4
				% DE/current-to-rand/1
				U(:, i) = V(:, i);
			else
				fprintf('BUG\n');
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
				k = selStrategy(i);				
				if fui < f(i)
					f(i) = fui;
					X(:, i) = U(:, i);
					CRMemory{k}(end + 1) = CR(k, i);
					if g <= LP
						ns(k, g) = ns(k, g) + 1;
					else
						ns(k, end) = ns(k, end) + 1;
					end
				else
					if g <= LP
						nf(k, g) = nf(k, g) + 1;
					else
						nf(k, end) = nf(k, end) + 1;
					end
				end
			end
		else
			for i = 1 : NP
				k = selStrategy(i);
				fui = feval(fitfun, U(:, i));
				counteval = counteval + 1;
				
				if fui < f(i)
					f(i) = fui;
					X(:, i) = U(:, i);
					CRMemory{k}(end + 1) = CR(k, i);
					if g <= LP
						ns(k, g) = ns(k, g) + 1;
					else
						ns(k, end) = ns(k, end) + 1;
					end
				else
					if g <= LP
						nf(k, g) = nf(k, g) + 1;
					else
						nf(k, end) = nf(k, end) + 1;
					end
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