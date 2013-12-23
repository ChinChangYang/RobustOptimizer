function [xmin, fmin, out] = dfdeeig(fitfun, lb, ub, maxfunevals, options)
% DFDEEIG "Eig" variant of Drift-free Differential Evolution
% DFDEEIG(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DFDEEIG(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 2 * numel(lb);
defaultOptions.maxfunevalsFactor =  0;
defaultOptions.R = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.F = 0.7;
defaultOptions.pmu = 0.5;
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
CR = options.CR;
F = options.F;
pmu = options.pmu;
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
		for i = 1 : NP
			r1 = floor(1 + NP * rand);
			r2 = floor(1 + NP * rand);
			
			while r1 == i
				r1 = floor(1 + NP * rand);
			end
			
			while r2 == i || r2 == r1
				r2 = floor(1 + NP * rand);
			end
			
			if rand <= pmu
				V(:, i) = X(:, i) + F * (X(:, r1) - X(:, r2));
			else
				sum2 = 0;
				while sum2 == 0
					r3 = floor(1 + NP * rand);
					r4 = floor(1 + NP * rand);
					
					while r3 == i
						r3 = floor(1 + NP * rand);
					end
					
					while r4 == i || r4 == r3
						r4 = floor(1 + NP * rand);
					end
					
					sum1 = 0;
					sum2 = 0;
					
					for j = 1 : D
						d1 = X(j, r1) - X(j, r2);
						d2 = X(j, r3) - X(j, r4) - 2 * X(j, i);
						sum1 = sum1 + d1 * d2;
						sum2 = sum2 + d2 * d2;
					end
				end
				
				K = sqrt(D) * sum1 / sum2;
				
				V(:, i) = X(:, i) + K * (X(:, r3) + X(:, r4) - 2 * X(:, i));
			end
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
				end
			end
		else
			counteval = counteval + NP;
			for i = 1 : NP
				fui = feval(fitfun, U(:, i));
				
				if fui < f(i)
					f(i) = fui;
					X(:, i) = U(:, i);
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