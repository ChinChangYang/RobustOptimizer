function [xmin, fmin, out] = cmaes(fitfun, lb, ub, maxfunevals, options)
% CMAES CMA-ES
% CMAES(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% CMAES(..., options) minimize the function by solver options.
%
% CMA-ES code is a modified version of that presented on
% (https://www.lri.fr/~hansen/cmaes_inmatlab.html#matlab).
if nargin <= 4
	options = [];
end

defaultOptions.NP = 100;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolStagnationIteration = Inf;
defaultOptions.usefunevals = inf;
defaultOptions.initial.xmin = [];
defaultOptions.initial.fmin = [];
defaultOptions.initial.m = [];
defaultOptions.initial.C = [];
defaultOptions.initial.ps = [];
defaultOptions.initial.pc = [];
defaultOptions.initial.sigma = [];
defaultOptions.initial.counteval = [];
defaultOptions.initial.countiter = [];
defaultOptions.ConstraintHandling = 'OnBound';
defaultOptions.EarlyStop = 'none';

options = setdefoptions(options, defaultOptions);
NP = options.NP;
usefunevals = options.usefunevals;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;

if isequal(options.ConstraintHandling, 'OnBound')
	onbound = true;
else
	onbound = false;
end

if ~isempty(strfind(options.EarlyStop, 'fitness'))
	EarlyStopOnFitness = true;
	AutoEarlyStop = false;
elseif ~isempty(strfind(options.EarlyStop, 'auto'))
	EarlyStopOnFitness = false;
	AutoEarlyStop = true;
else
	EarlyStopOnFitness = false;
	AutoEarlyStop = false;
end

if ~isempty(options.initial)
	options.initial = setdefoptions(options.initial, defaultOptions.initial);
	xmin		= options.initial.xmin;
	fmin		= options.initial.fmin;
	m			= options.initial.m;
	C			= options.initial.C;
	ps			= options.initial.ps;
	pc			= options.initial.pc;
	sigma		= options.initial.sigma;
	counteval	= options.initial.counteval;
	countiter	= options.initial.countiter;
else
	xmin		= [];
	fmin		= [];
	m			= [];
	C			= [];
	ps			= [];
	pc			= [];
	sigma		= [];
	counteval	= [];
	countiter	= [];
end

D = numel(lb);

% Initialize variables
out = initoutput(RecordPoint, D, NP, maxfunevals);

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
end

% counteval
if isempty(counteval)
	counteval = 0;
end

% countiter
if isempty(countiter)
	countiter = 1;
end

% fmin
if isempty(fmin)
	fmin = inf;
end

% m
if isempty(m)
	m = (lb + ub) / 2;
end

% C
if isempty(C)
	diagD = ones(D, 1);
	C = diag(diagD.^2);
end

% ps
if isempty(ps)
	ps = zeros(D, 1);
end

% pc
if isempty(pc)
	pc = zeros(D, 1);
end

% sigma
if isempty(sigma)	
	sigma = max(ub - lb) / 8;
end

% Covariance matrix repair
C = triu(C) + triu(C,1)';
C = real(C);
[B, temp] = eig(C);
B = real(B);
temp = real(temp);
diagD = diag(temp);

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

% Initialize variables
BD = B .* repmat(diagD', D, 1);
mu = floor(0.5 * NP);
w = log(mu + 0.5)-log(1:mu)';
mueff = sum(w)^2 / sum(w.^2);
w = w / sum(w);
cs = (mueff + 2) / (D + mueff + 3);
ds = 1 + 2 * max(0, sqrt((mueff - 1) / (D + 1)) - 1) + cs;
cc = (4 + mueff / D) / (D + 4 + 2 * mueff / D);
c1 = 2 / ((D + 1.3)^2 + mueff);
cmu = 2 * (mueff - 2 + 1 / mueff) / ((D + 2)^2 + mueff);
chiN = D ^ 0.5 * (1 - 1 / (4 * D) + 1 / (21 * D ^ 2));

while true
	% Generate population
	Z = randn(D, NP);
	X = repmat(m, 1, NP) + sigma * BD * Z;
	XVALID = X;
	
	% Correction for outside of boundaries	
	if onbound
		for i = 1 : NP
			for j = 1 : D
				if XVALID(j, i) < lb(j)
					XVALID(j, i) = lb(j);
				elseif XVALID(j, i) > ub(j)
					XVALID(j, i) = ub(j);
				end
			end
		end
	end
	
	% Evaluation
	fx = zeros(1, NP);
	for i = 1 : NP
		fx(i) = feval(fitfun, XVALID(:, i));
		counteval = counteval + 1;
	end
	
	for i = 1 : NP		
		penalty = 1e10 * std(fx) * norm(XVALID(:, i) - X(:, i));
		fx(i) = fx(i) + penalty;
	end
	
	% Sort	
	[fx, fidx] = sort(fx);
	X = X(:, fidx);
	Z = Z(:, fidx);
	
	if fx(1) < fmin
		fmin = fx(1);
		xmin = X(:, 1);
	end
	
	% Update
	zmean = Z(:, 1 : mu) * w;
	mold = m;
	m = X(:, 1 : mu) * w;
	ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * (B * zmean);
	hsig = norm(ps) ...
		/ sqrt(1 - (1 - cs)^(2 * countiter)) ...
		/ chiN < 1.4 + 2 / (D + 1);
	pc = (1 - cc) * pc ...
		+ hsig * sqrt(cc * (2 - cc) * mueff) * (m - mold) / sigma;
	arpos = (X(:, 1 : mu) - repmat(mold, 1, mu)) / sigma;
	C = (1 - c1 - cmu + (1 - hsig) * c1 * cc * (2 - cc)) * C ...
		+ c1 * (pc * pc') ...
		+ cmu ...
		* arpos * (repmat(w, 1, D) .* arpos');
	sigma = sigma * exp(min(1, (sqrt(sum(ps.^2)) / chiN - 1) * cs / ds));
	
	% Covariance matrix repair
	C = triu(C) + triu(C,1)';
	C = real(C);
	[B, temp] = eig(C);
	B = real(B);
	temp = real(temp);
	diagD = diag(temp);
	
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
	
	% Correction for outside of boundaries	
	if onbound
		for j = 1 : D
			if m(j) < lb(j)
				m(j) = lb(j);
			elseif m(j) > ub(j)
				m(j) = ub(j);
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			X, X, fx, countiter, XX, YY, ZZ);
	end
	
	% Record
	out = updateoutput(out, X, fx, counteval, countiter);
	
	% Iteration counter
	countiter = countiter + 1;

	% Termination conditions
	outofmaxfunevals = counteval > maxfunevals - NP;
	outofusefunevals = counteval > usefunevals - NP;
	if ~EarlyStopOnFitness && ~AutoEarlyStop
		if outofmaxfunevals || outofusefunevals
			break;
		end
	elseif AutoEarlyStop
		reachftarget = min(fx) <= ftarget;
		TolX = 10 * eps(mean(X(:)));
		solutionconvergence = std(X(:)) <= TolX;
		TolFun = 10 * eps(mean(fx));
		functionvalueconvergence = std(fx(:)) <= TolFun;
		
		if outofmaxfunevals || ...
				reachftarget || ...
				solutionconvergence || ...
				functionvalueconvergence
			break;
		end
	elseif EarlyStopOnFitness
		reachftarget = min(fx) <= ftarget;
		
		if outofmaxfunevals || ...
				reachftarget
			break;
		end
	end
end

final.xmin = xmin;
final.fmin = fmin;
final.m = m;
final.C = C;
final.ps = ps;
final.pc = pc;
final.sigma = sigma;
final.counteval = counteval;
final.countiter = countiter;

out = finishoutput(out, X, fx, counteval, countiter, ...
	'final', final);
end
