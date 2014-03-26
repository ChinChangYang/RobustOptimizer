function [xmin, fmin, out] = jadeeig(fitfun, lb, ub, maxfunevals, options)
% JADEEIG JADE/eig
% JADEEIG(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% JADEEIG(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.R = 0.5;
defaultOptions.mu_CR = 0.5;
defaultOptions.mu_F = 0.7;
defaultOptions.delta_CR = 0.1;
defaultOptions.delta_F = 0.1;
defaultOptions.p = 0.05;
defaultOptions.w = 0.1;
defaultOptions.Display = 'off';
defaultOptions.FactorNP = 1;
defaultOptions.Restart = 0;
defaultOptions.RecordPoint = 100;
defaultOptions.Noise = false;
defaultOptions.SampleFactor = 1.01;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.TolStagnationIteration = 30;
defaultOptions.X = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = max(1, options.dimensionFactor);
R = options.R;
mu_CR = options.mu_CR;
mu_F = options.mu_F;
delta_CR = options.delta_CR;
delta_F = options.delta_F;
p = options.p;
w = options.w;
isDisplayIter = strcmp(options.Display, 'iter');
FactorNP = max(1, options.FactorNP);
Restart = max(0, round(options.Restart));
RecordPoint = max(1, floor(options.RecordPoint));
Noise = options.Noise;
SampleFactor = options.SampleFactor;
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
TolStagnationIteration = options.TolStagnationIteration;
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
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
	'mu_F', 'mu_CR');
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
	pbest_size = p * NP;	
	countStagnation = 0;
	
	% Evaluation
	f = zeros(1, NP);
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
	
	[f, fidx] = sort(f);
	X = X(:, fidx);
	
	% Initialize archive
	A = zeros(D, 2 * NP);
	A(:, 1 : NP) = X;
	
	% Display
	if isDisplayIter
		displayitermessages(...
			X, U, f, countiter, XX, YY, ZZ, 'mu_F', mu_F, 'mu_CR', mu_CR);
	end
	
	% Record
	out = updateoutput(out, X, f, counteval, ...
		'mu_F', mu_F, 'mu_CR', mu_CR);
	
	% Iteration counter
	countiter = countiter + 1;
	
	while true
		% Termination conditions
		outofmaxfunevals = counteval > maxfunevals - NP;
		reachftarget = min(f) <= ftarget;
		fitnessconvergence = isConverged(f, TolFun);
		solutionconvergence = isConverged(X, TolX);
		stagnation = countStagnation >= TolStagnationIteration;
		
		% Convergence conditions
		if outofmaxfunevals || reachftarget || fitnessconvergence || ...
				solutionconvergence || stagnation
			break;
		end
		
		% Scaling factor and crossover rate
		S_F = zeros(1, NP);
		S_CR = zeros(1, NP);
		CR = mu_CR + delta_CR * randn(1, NP);
		CR(CR > 1) = 1;
		CR(CR < 0) = 0;
		F = cauchyrnd(mu_F, delta_F, NP, 1);
		F(F > 1) = 1;
		
		for retry = 1 : 3
			if all(F > 0)
				break;
			end
			
			F(F <= 0) = cauchyrnd(mu_F, delta_F, sum(F <= 0), 1);
			F(F > 1) = 1;
		end
		
		F(F <= 0) = eps;
		
		A_Counter = 0;
		XA = [X, A];
		
		% Mutation
		for i = 1 : NP
			
			for checkRetry = 1 : 3
				
				% Generate pbest_idx
				for retry = 1 : 3
					pbest_idx = max(1, ceil(rand * pbest_size));
					if ~all(X(:, pbest_idx) == X(:, i))
						break;
					end
				end
				
				% Generate r1
				for retry = 1 : 3
					r1 = floor(1 + NP * rand);
					if i ~= r1
						break;
					end
				end
				
				% Generate r2
				for retry = 1 : 3
					r2 = floor(1 + 2 * NP * rand);
					if ~(all(X(:, i) == XA(:, r2)) || all(X(:, r1) == XA(:, r2)))
						break;
					end
				end
				
				V(:, i) = X(:, i) + F(i) * (X(:, pbest_idx) - X(:, i) + X(:, r1) - XA(:, r2));
				
				% Check boundary
				if all(V(:, i) > lb) && all(V(:, i) < ub)
					break;
				end
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
					if rand < CR(i) || j == jrand
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
					if rand < CR(i) || j == jrand
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
			displayitermessages(...
				X, U, f, countiter, XX, YY, ZZ, 'mu_F', mu_F, 'mu_CR', mu_CR);
		end
		
		% Selection
		FailedIteration = true;
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
					A(:, NP + A_Counter + 1) = U(:, i);
					S_CR(A_Counter + 1) = CR(i);
					S_F(A_Counter + 1) = F(i);
					A_Counter = A_Counter + 1;
					FailedIteration = false;
				end
			end
		else
			for i = 1 : NP
				fui = feval(fitfun, U(:, i));
				counteval = counteval + 1;
				
				if fui < f(i)
					A(:, NP + A_Counter + 1) = X(:, i);
					f(i) = fui;
					X(:, i) = U(:, i);
					S_CR(A_Counter + 1) = CR(i);
					S_F(A_Counter + 1) = F(i);
					A_Counter = A_Counter + 1;
					FailedIteration = false;
				end
			end
		end
		
		% Update archive
		rand_idx = randperm(NP + A_Counter);
		A(:, 1 : NP) = A(:, rand_idx(1 : NP));
		
		% Update CR and F
		if A_Counter > 0
			mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : A_Counter));
			mu_F = (1-w) * mu_F + w * sum(S_F(1 : A_Counter).^2) / sum(S_F(1 : A_Counter));
		end
		
		% Sort
		[f, fidx] = sort(f);
		X = X(:, fidx);
		
		% Record
		out = updateoutput(out, X, f, counteval, ...
			'mu_F', mu_F, 'mu_CR', mu_CR);
		
		% Iteration counter
		countiter = countiter + 1;
		
		% Stagnation iteration
		if FailedIteration
			countStagnation = countStagnation + 1;
		else
			countStagnation = 0;
		end
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

out = finishoutput(out, X, f, counteval, ...
	'mu_F', mu_F, 'mu_CR', mu_CR);
fmin = out.bestever.fmin;
xmin = out.bestever.xmin;
end
