function [xmin, fmin, out] = jadeeigv2(fitfun, lb, ub, maxfunevals, options)
% JADEEIGV2 JADE/eig with F<0 exception handling by F=mu_F
% JADEEIGV2(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% JADEEIGV2(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.maxfunevalsFactor = 0;
defaultOptions.R = 0.7;
defaultOptions.mu_CR = 0.5;
defaultOptions.mu_F = 0.5;
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

options = setdefoptions(options, defaultOptions);
dimensionFactor = max(1, options.dimensionFactor);
maxfunevalsFactor = options.maxfunevalsFactor;
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

D = numel(lb);
NP = max(ceil(dimensionFactor * D), floor(maxfunevalsFactor * maxfunevals));

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
	
	% Record minimal function values
	out.fmin = inf(1, RecordPoint);
	out.fmin(end) = min(f);
	
	% Record mean function values
	out.fmean = inf(1, RecordPoint);
	out.fmean(end) = mean(f);
	
	% Record std. function values
	out.fstd = inf(1, RecordPoint);
	out.fstd(end) = std(f);
	
	% Record minimizer
	out.xmin = inf(D, RecordPoint);
	out.xmin(:, end) = xmin;
	
	% Record mean of minimizers
	out.xmean = inf(D, RecordPoint);
	out.xmean(:, end) = mean(X, 2);
	
	% Record std. of minimizers
	out.xstd = inf(D, RecordPoint);
	out.xstd(:, end) = std(X, 0, 2);
	
	% Record function evaluations
	out.fes = zeros(1, RecordPoint);
	out.fes(end) = maxfunevals;
	
	% Record mean of samples
	out.samplemean = zeros(1, RecordPoint);
	out.samplemean(end) = 1;
	
	% Record population distances
	distance = sampledistance(X);
	out.distancemin = zeros(1, RecordPoint);
	out.distancemin(end) = min(distance);
	out.distancemax = zeros(1, RecordPoint);
	out.distancemax(end) = max(distance);
	out.distancemean = zeros(1, RecordPoint);
	out.distancemean(end) = mean(distance);
	out.distancemedian = zeros(1, RecordPoint);
	out.distancemedian(end) = median(distance);
	
	% Record scale factors
	out.mu_F = zeros(1, RecordPoint);
	out.mu_F(end) = mu_F;
	
	return;
end

% Initialize variables
counteval = 0;
countiter = 1;
recordFEs = linspace(0, 1, RecordPoint) * maxfunevals;
out.fmin = inf(1, RecordPoint);
out.fmean = inf(1, RecordPoint);
out.fstd = inf(1, RecordPoint);
out.xmean = inf(D, RecordPoint);
out.xmin = inf(D, RecordPoint);
out.xstd = inf(D, RecordPoint);
out.fes = zeros(1, RecordPoint);
out.samplemean = zeros(1, RecordPoint);
out.distancemin = zeros(1, RecordPoint);
out.distancemax = zeros(1, RecordPoint);
out.distancemean = zeros(1, RecordPoint);
out.distancemedian = zeros(1, RecordPoint);
out.mu_F = zeros(1, RecordPoint);
out.bestever.fmin = Inf;
iRecordFEs = 1;
initNP = NP;
MAX_NP = 2000;

% Initialize contour data
if isDisplayIter
	vx = linspace(lb(1), ub(1), 50);
	vy = linspace(lb(2), ub(2), 50);
	[XX, YY] = meshgrid(vx, vy);
	ZZ = zeros(numel(vx), numel(vy));
	for i = 1 : numel(vx)
		for j = 1 : numel(vy)
			ZZ(i, j) = feval(fitfun, [XX(i, j); YY(i, j)]);
		end
	end
	ZZ = log(ZZ - min(ZZ(:)) + 1);
end

for iRestart = 1 : (Restart + 1)
	NP = min(MAX_NP, round(initNP * FactorNP ^ (iRestart - 1)));
	
	% Initialize population
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
	
	% Initialize variables
	V = X;
	U = X;
	XT = X;
	VT = X;
	UT = X;
	realSamples = 1;
	pbest_size = p * NP;	
	
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
		displayitermessages(X, U, f, countiter, XX, YY, ZZ);
	end
	
	% Record minimal function values
	while counteval >= recordFEs(iRecordFEs)
		[fmin, fminidx] = min(f);
		xmin = X(:, fminidx);
		distance = sampledistance(X);
		out.fmin(iRecordFEs) = fmin;
		out.fmean(iRecordFEs) = mean(f);
		out.fstd(iRecordFEs) = std(f);
		out.xmin(:, iRecordFEs) = xmin;
		out.xmean(:, iRecordFEs) = mean(X, 2);
		out.xstd(:, iRecordFEs) = std(X, 0, 2);
		out.fes(iRecordFEs) = counteval;
		out.samplemean(iRecordFEs) = realSamples;
		out.distancemin(iRecordFEs) = min(distance);
		out.distancemax(iRecordFEs) = max(distance);
		out.distancemean(iRecordFEs) = mean(distance);
		out.distancemedian(iRecordFEs) = median(distance);
		out.mu_F(iRecordFEs) = mu_F;
		iRecordFEs = iRecordFEs + 1;
		if iRecordFEs > RecordPoint
			break;
		end
	end
	
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
		
		F(F <= 0) = mu_F;
		
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
					U(j, i) = 2 * lb(j) - U(j, i);
				elseif U(j, i) > ub(j)
					U(j, i) = 2 * ub(j) - U(j, i);
				else
					continue;
				end
				
				if U(j, i) < lb(j)
					U(j, i) = lb(j);
				elseif U(j, i) > ub(j)
					U(j, i) = ub(j);
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
					A(:, NP + A_Counter + 1) = U(:, i);
					S_CR(A_Counter + 1) = CR(i);
					S_F(A_Counter + 1) = F(i);
					A_Counter = A_Counter + 1;
				end
			end
		else
			for i = 1 : NP
				fui = feval(fitfun, U(:, i));
				counteval = counteval + 1;
				
				if fui < f(i)
					f(i) = fui;
					X(:, i) = U(:, i);
					A(:, NP + A_Counter + 1) = U(:, i);
					S_CR(A_Counter + 1) = CR(i);
					S_F(A_Counter + 1) = F(i);
					A_Counter = A_Counter + 1;
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
		while counteval >= recordFEs(iRecordFEs)
			[fmin, fminidx] = min(f);
			xmin = X(:, fminidx);
			distance = sampledistance(X);
			out.fmin(iRecordFEs) = fmin;
			out.fmean(iRecordFEs) = mean(f);
			out.fstd(iRecordFEs) = std(f);
			out.xmin(:, iRecordFEs) = xmin;
			out.xmean(:, iRecordFEs) = mean(X, 2);
			out.xstd(:, iRecordFEs) = std(X, 0, 2);
			out.fes(iRecordFEs) = counteval;
			out.samplemean(iRecordFEs) = realSamples;
			out.distancemin(iRecordFEs) = min(distance);
			out.distancemax(iRecordFEs) = max(distance);
			out.distancemean(iRecordFEs) = mean(distance);
			out.distancemedian(iRecordFEs) = median(distance);
			out.mu_F(iRecordFEs) = mu_F;
			iRecordFEs = iRecordFEs + 1;
			if iRecordFEs > RecordPoint
				break;
			end
		end
		
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

distance = sampledistance(X);
out.fmin(iRecordFEs:end) = min(f);
out.fmean(iRecordFEs:end) = mean(f);
out.fstd(iRecordFEs:end) = std(f);
nRemaining = RecordPoint - iRecordFEs + 1;
out.xmin(:, iRecordFEs:RecordPoint) = repmat(xmin, 1, nRemaining);
out.xmean(:, iRecordFEs:RecordPoint) = repmat(mean(X, 2), 1, nRemaining);
out.xstd(:, iRecordFEs:RecordPoint) = repmat(std(X, 0, 2), 1, nRemaining);
out.fes(iRecordFEs:RecordPoint) = counteval;
out.samplemean(iRecordFEs:RecordPoint) = realSamples;
out.distancemin(iRecordFEs:RecordPoint) = min(distance);
out.distancemax(iRecordFEs:RecordPoint) = max(distance);
out.distancemean(iRecordFEs:RecordPoint) = mean(distance);
out.distancemedian(iRecordFEs:RecordPoint) = median(distance);
out.mu_F(iRecordFEs:RecordPoint) = mu_F;

fmin = out.bestever.fmin;
xmin = out.bestever.xmin;
end
