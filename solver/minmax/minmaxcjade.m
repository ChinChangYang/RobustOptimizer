function [xminmax1, xminmax2, fminmax, out] = minmaxcjade(fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXCJADE Coevolutionary JADE for min-max problems
% MINMAXCJADE(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MINMAXCJADE(..., options1) minimizes the function with the given options
% Options1 for the 1st layer.
% MINMAXCJADE(..., options1, options2) minimize the function with the given
% options Options2 for the 2nd layer.
if nargin <= 6
	options1 = [];
end

if nargin <= 7
	options2 = [];
end

D1 = numel(lb1);
D2 = numel(lb2);

% Default options for Layer 1
defaultOptions1.dimensionFactor = 5;
defaultOptions1.mu_CR = 0.5;
defaultOptions1.mu_F = 1.0;
defaultOptions1.delta_CR = 0.1;
defaultOptions1.delta_F = 0.1;
defaultOptions1.p = 0.05;
defaultOptions1.w = 0.1;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.ftarget = -Inf;
defaultOptions1.TolFun = eps;
defaultOptions1.TolX = 100 * eps;
options1 = setdefoptions(options1, defaultOptions1);
NP1 = ceil(options1.dimensionFactor * D1);
RecordPoint1 = options1.RecordPoint;
isDisplayIter = strcmp(options1.Display, 'iter');
TolFun1 = options1.TolFun;

% Default options for Layer 2
defaultOptions2 = defaultOptions1;
defaultOptions2.RecordPoint = 0;
options2 = setdefoptions(options2, defaultOptions2);
NP2 = ceil(options2.dimensionFactor * D2);
TolFun2 = options2.TolFun;

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint1, D1, NP1, maxfunevals, 'mu_F', 'mu_CR');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = minmaxcontourdata(D1, lb1, ub1, lb2, ub2, fitfun);
end

% Initialize population 1
if NP1 < 1e1
	LHS1 = lhsdesign(NP1, D1, 'iteration', 10)';
elseif NP1 < 1e2
	LHS1 = lhsdesign(NP1, D1, 'iteration', 2)';
else
	LHS1 = rand(D1, NP1);
end

X1 = zeros(D1, NP1);
for i = 1 : NP1
	X1(:, i) = lb1 + (ub1 - lb1) .* LHS1(:, i);
end

V1 = X1;
U1 = X1;
A1 = zeros(D1, 2 * NP1);
A1(:, 1 : NP1) = X1;
mu_F1 = options1.mu_F;
delta_F1 = options1.delta_F;
mu_CR1 = options1.mu_CR;
delta_CR1 = options1.delta_CR;
pbest1 = options1.p * NP1;

% Initialize population 2
if NP2 < 1e1
	LHS2 = lhsdesign(NP2, D2, 'iteration', 10)';
elseif NP1 < 1e2
	LHS2 = lhsdesign(NP2, D2, 'iteration', 2)';
else
	LHS2 = rand(D2, NP2);
end

X2 = zeros(D2, NP2);
for i = 1 : NP2
	X2(:, i) = lb2 + (ub2 - lb2) .* LHS2(:, i);
end

V2 = X2;
U2 = X2;
A2 = zeros(D2, 2 * NP2);
A2(:, 1 : NP2) = X2;
mu_F2 = options2.mu_F;
delta_F2 = options2.delta_F;
mu_CR2 = options2.mu_CR;
delta_CR2 = options2.delta_CR;
pbest2 = options2.p * NP2;

% Evaluation
f = zeros(NP1, NP2);
for i = 1 : NP1
	for j = 1 : NP2
		f(i, j) = feval(fitfun, X1(:, i), X2(:, j));
		counteval = counteval + 1;
	end
end

fmax = max(f, [], 2);
[fmax, fmaxidx] = sort(fmax);
X1 = X1(:, fmaxidx);

fmin = min(f);
[fmin, fminidx] = sort(fmin, 'descend');
X2 = X2(:, fminidx);

f = f(fmaxidx, fminidx);

% Display
if isDisplayIter	
	S = zeros(D1 + D2, NP1 * NP2);
	for i = 1 : NP1
		for j = 1 : NP2
			S(:, (i - 1)* NP2 + j) = [X1(:, i); X2(:, j)];
		end
	end
	displayitermessages(S, S, fmax, countiter, XX, YY, ZZ, 'mu_F1', mu_F1, 'mu_CR1', mu_CR1);
end

out = updateoutput(out, X1, fmax, counteval, 'mu_F', mu_F1, 'mu_CR', mu_CR1);
countiter = countiter + 1;

while true
	% Termination conditions
	stdfmax = std(fmax);
	stdfmin = std(fmin);
	stdX1 = std(X1, 0, 2);
	meanX1 = mean(X1, 2);
	stdX2 = std(X2, 0, 2);
	meanX2 = mean(X2, 2);
	outofmaxfunevals = counteval > maxfunevals - NP1 * NP2;
	fitnessconvergence1 = stdfmax <= mean(fmax) * 100 * eps || stdfmax <= realmin || stdfmax <= TolFun1;
	fitnessconvergence2 = stdfmin <= mean(fmin) * 100 * eps || stdfmin <= realmin || stdfmin <= TolFun2;
	solutionconvergence1 = all(stdX1 <= meanX1 * 100 * eps) || all(stdX1 <= 100 * realmin) || ...
		all(stdX1 <= options1.TolX);
	solutionconvergence2 = all(stdX2 <= meanX2 * 100 * eps) || all(stdX2 <= 100 * realmin) || ...
		all(stdX2 <= options2.TolX);
	
	if outofmaxfunevals || fitnessconvergence1 || fitnessconvergence2 || ...
			solutionconvergence1 || solutionconvergence2
		break;
	end
	
	% X1's Scaling factor and crossover rate
	S_F1 = zeros(1, NP1);
	S_CR1 = zeros(1, NP1);
	CR1 = mu_CR1 + delta_CR1 * randn(1, NP1);
	CR1(CR1 > 1) = 1;
	CR1(CR1 < 0) = 0;
	F1 = cauchyrnd(mu_F1, delta_F1, NP1, 1);
	F1(F1 > 1) = 1;
	
	for retry = 1 : 3
		if all(F1 > 0)
			break;
		end
		
		F1(F1 <= 0) = cauchyrnd(mu_F1, delta_F1, sum(F1 <= 0), 1);
		F1(F1 > 1) = 1;
	end
	
	F1(F1 <= 0) = 100 * eps * (1 + rand);
		
	A1_Counter = 0;
	XA1 = [X1, A1];
	
	% X1's Mutation
	for i = 1 : NP1
		for checkRetry = 1 : 3
			% Generate pbest_idx
			for retry = 1 : 3
				pbest1_idx = max(1, ceil(rand * pbest1));
				if ~all(X1(:, pbest1_idx) == X1(:, i))
					break;
				end
			end
			
			% Generate r1
			for retry = 1 : 3
				r1 = floor(1 + NP1 * rand);
				if i ~= r1
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : 3
				r2 = floor(1 + 2 * NP1 * rand);
				if ~(all(X1(:, i) == XA1(:, r2)) || all(X1(:, r1) == XA1(:, r2)))
					break;
				end
			end
			
			V1(:, i) = X1(:, i) + F1(i) * (X1(:, pbest1_idx) - X1(:, i) + X1(:, r1) - XA1(:, r2));
			
			% Check boundary
			if all(V1(:, i) > lb1) && all(V1(:, i) < ub1)
				break;
			end
		end
	end
	
	for i = 1 : NP1
		% Binominal Crossover
		jrand = floor(1 + D1 * rand);
		for j = 1 : D1
			if rand < CR1(i) || j == jrand
				U1(j, i) = V1(j, i);
			else
				U1(j, i) = X1(j, i);
			end
		end
	end
	
	% X1's Repair
	for i = 1 : NP1
		for j = 1 : D1
			if U1(j, i) < lb1(j)
				U1(j, i) = 2 * lb1(j) - U1(j, i);
			elseif U1(j, i) > ub1(j)
				U1(j, i) = 2 * ub1(j) - U1(j, i);
			else
				continue;
			end
			
			if U1(j, i) < lb1(j)
				U1(j, i) = lb1(j);
			elseif U1(j, i) > ub1(j)
				U1(j, i) = ub1(j);
			end
		end
	end
	
	% X2's Scaling factor and crossover rate
	S_F2 = zeros(1, NP2);
	S_CR2 = zeros(1, NP2);
	CR2 = mu_CR2 + delta_CR2 * randn(1, NP2);
	CR2(CR2 > 1) = 1;
	CR2(CR2 < 0) = 0;
	F2 = cauchyrnd(mu_F2, delta_F2, NP2, 1);
	F2(F2 > 1) = 1;
	
	for retry = 1 : 3
		if all(F2 > 0)
			break;
		end
		
		F2(F2 <= 0) = cauchyrnd(mu_F2, delta_F2, sum(F2 <= 0), 1);
		F2(F2 > 1) = 1;
	end
	
	F2(F2 <= 0) = 100 * eps * (1 + rand);
		
	A2_Counter = 0;
	XA2 = [X2, A2];
	
	% X1's Mutation
	for i = 1 : NP2
		for checkRetry = 1 : 3
			% Generate pbest_idx
			for retry = 1 : 3
				pbest2_idx = max(1, ceil(rand * pbest2));
				if ~all(X2(:, pbest2_idx) == X2(:, i))
					break;
				end
			end
			
			% Generate r1
			for retry = 1 : 3
				r1 = floor(1 + NP2 * rand);
				if i ~= r1
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : 3
				r2 = floor(1 + 2 * NP2 * rand);
				if ~(all(X2(:, i) == XA2(:, r2)) || all(X2(:, r1) == XA2(:, r2)))
					break;
				end
			end
			
			V2(:, i) = X2(:, i) + F2(i) * (X2(:, pbest2_idx) - X2(:, i) + X2(:, r1) - XA2(:, r2));
			
			% Check boundary
			if all(V2(:, i) > lb2) && all(V2(:, i) < ub2)
				break;
			end
		end
	end
	
	for i = 1 : NP2
		% Binominal Crossover
		jrand = floor(1 + D2 * rand);
		for j = 1 : D2
			if rand < CR2(i) || j == jrand
				U2(j, i) = V2(j, i);
			else
				U2(j, i) = X2(j, i);
			end
		end
	end
	
	% X2's Repair
	for i = 1 : NP2
		for j = 1 : D2
			if U2(j, i) < lb2(j)
				U2(j, i) = 2 * lb2(j) - U2(j, i);
			elseif U2(j, i) > ub2(j)
				U2(j, i) = 2 * ub2(j) - U2(j, i);
			else
				continue;
			end
			
			if U2(j, i) < lb2(j)
				U2(j, i) = lb2(j);
			elseif U2(j, i) > ub2(j)
				U2(j, i) = ub2(j);
			end
		end
	end
	
	% Evaluation
	fu = zeros(NP1, NP2);
	for i = 1 : NP1
		for j = 1 : NP2
			fu(i, j) = feval(fitfun, U1(:, i), U2(:, j));
			counteval = counteval + 1;
		end
	end
	
	fumax = max(fu, [], 2);
	fumin = min(fu);
	
	% X1's Selection
	for i = 1 : NP1
		if fumax(i) < fmax(i)
			fmax(i) = fumax(i);
			X1(:, i) = U1(:, i);
			A1(:, NP1 + A1_Counter) = U1(:, i);
			S_CR1(A1_Counter + 1) = CR1(i);
			S_F1(A1_Counter + 1) = F1(i);
			A1_Counter = A1_Counter + 1;
		end
	end
	
	% X1's archive update
	rand_idx = randperm(NP1 + A1_Counter);
	A1(:, 1 : NP1) = A1(:, rand_idx(1 : NP1));
	
	% X1's CR and F update
	if A1_Counter > 0
		mu_CR1 = (1-options1.w) * mu_CR1 + options1.w * mean(S_CR1(1 : A1_Counter));
		mu_F1 = (1-options1.w) * mu_F1 + options1.w * sum(S_F1(1 : A1_Counter).^2) / sum(S_F1(1 : A1_Counter));
	end
	
	% X2's Selection
	for i = 1 : NP2
		if fumin(i) > fmin(i)
			fmin(i) = fumin(i);
			X2(:, i) = U2(:, i);
			A2(:, NP2 + A2_Counter) = U2(:, i);
			S_CR2(A2_Counter + 1) = CR2(i);
			S_F2(A2_Counter + 1) = F2(i);
			A2_Counter = A2_Counter + 1;
		end
	end
	
	% X2's archive update
	rand_idx = randperm(NP2 + A2_Counter);
	A2(:, 1 : NP2) = A2(:, rand_idx(1 : NP2));
	
	% X2's CR and F update
	if A2_Counter > 0
		mu_CR2 = (1-options2.w) * mu_CR2 + options2.w * mean(S_CR2(1 : A2_Counter));
		mu_F2 = (1-options2.w) * mu_F2 + options2.w * sum(S_F2(1 : A2_Counter).^2) / sum(S_F2(1 : A2_Counter));
	end
	
	% Display
	if isDisplayIter
		S = zeros(D1 + D2, NP1 * NP2);
		for i = 1 : NP1
			for j = 1 : NP2
				S(:, (i - 1)* NP2 + j) = [X1(:, i); X2(:, j)];
			end
		end
		displayitermessages(S, S, fmax, countiter, XX, YY, ZZ, 'mu_F1', mu_F1, 'mu_CR1', mu_CR1);
	end
		
	% Evaluation
	for i = 1 : NP1
		for j = 1 : NP2
			f(i, j) = feval(fitfun, X1(:, i), X2(:, j));
			counteval = counteval + 1;
		end
	end
	
	fmax = max(f, [], 2);
	[fmax, fmaxidx] = sort(fmax);
	X1 = X1(:, fmaxidx);
	
	fmin = min(f);
	[fmin, fminidx] = sort(fmin, 'descend');
	X2 = X2(:, fminidx);
	
	f = f(fmaxidx, fminidx);	
	
	out = updateoutput(out, X1, fmax, counteval, 'mu_F', mu_F1, 'mu_CR', mu_CR1);
	countiter = countiter + 1;
end

xminmax1 = X1(:, 1);
xminmax2 = X2(:, 1);
fminmax = fmax(1);
out = finishoutput(out, X1, fmax, counteval, 'mu_F', mu_F1, 'mu_CR', mu_CR1);
end