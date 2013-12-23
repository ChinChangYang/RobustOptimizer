function [xmin, ...      % minimum search point of last iteration
	fmin, ...      % function value of xmin
	counteval, ... % number of function evaluations done
	stopflag, ...  % stop criterion reached
	out, ...     % struct with various histories and solutions
	bestever ... % struct containing overall best solution (for convenience)
	] = kmeansde( ...
	fitfun, ...    % name of objective/fitness function
	xstart, ...    % objective variables initial point, determines N
	insteps, ...   % initial coordinate wise population steps
	inopts, ...    % options struct, see defopts below
	varargin )     % arguments passed to objective function
%KMEANSDE K-means Differential Evolution

% check input arguments
if nargin < 1
	% self-test
	clc;
	clear all;
	fitfun = 'fgeneric';
	xstart = zeros(2, 1);
	insteps = 10 * ones(2, 1);
	varargin = {};
	inopts.maxfunevals = 1e6;
	inopts.plotting = true;
	inopts.de_more_plot = true;
	
	s = RandStream('mt19937ar', 'Seed', 14);
	RandStream.setDefaultStream(s);
	fgeneric('finalize');
	datapath = 'out/bbobdata';
	Foptions.algName = 'kmeansde';
	countfitness = 0;
	countmaxfunevals = 0;
	counttolfun = 0;
	allinst = 4;
	
	for inst = allinst
		inopts.ftarget = fgeneric('initialize', 23, inst, datapath, Foptions);
		[xmin, fmin, counteval, stopflag, out, bestever] = ...
			kmeansde(fitfun, xstart, insteps, inopts, varargin{:});
		
		if strcmp('fitness', stopflag)
			countfitness = countfitness + 1;
		elseif strcmp('maxfunevals', stopflag)
			countmaxfunevals = countmaxfunevals + 1;
		elseif strcmp('TolFun', stopflag)
			counttolfun = counttolfun + 1;
		end
		
		fgeneric('finalize');
	end
	
	fprintf('success rate: %.0f%%\n', 100 * countfitness / length(allinst));
	fprintf('out of time: %.0f%%\n', 100 * countmaxfunevals / length(allinst));
	fprintf('TolFun rate: %.0f%%\n', 100 * counttolfun / length(allinst));
	
	return;
end

if nargin < 2
	xstart = [];
end

if isempty(xstart)
	error('Initial search point, and problem dimension, not determined');
end

if nargin < 3
	insteps = [];
end

if isa(insteps, 'struct') || isempty(insteps)
	error(['Third argument INSTEPS must be (or eval to) a scalar ' ...
		'or a column vector of size(X0,1)']);
end

if any(insteps == 0)
	error('some of insteps are zero.');
end

D = length(xstart);
defopts = desacsdefoptions(D);

if nargin < 4 || isempty(inopts)
	opts = defopts;
else
	opts = setdefoptions(inopts, defopts);
end

% initialize variables
counteval = 0;
countiter = 0;
stopflag = [];
NP = round(opts.de_globalpopsize);
ftarget = opts.ftarget;
% K = ceil(NP / D / 6);
K = min(2 * D, ceil(NP / D / 6));
F = opts.de_f;
L = opts.contourl;
R = opts.contourr;
U = opts.contouru;
D = length(xstart);
contourx = [];
contoury = [];
contourz = [];

% create population
x = repmat(xstart - 0.5 * insteps, 1, NP) + diag(insteps) * rand(D, NP);
f = zeros(1, NP);

% evaluate population
countiter = countiter + 1;
for i = 1 : NP
	f(i) = feval(fitfun, x(:, i));
	counteval = counteval + 1;
end

[f, idx] = sort(f);
x = x(:, idx);
xmin = x(:, 1);
fmin = f(1);
bestever.xmin = xmin;
bestever.fmin = fmin;

mindist = Inf;
maxdist = 0;
for i = 2 : NP
	dist = norm(xmin - x(:, i));
	if dist < mindist
		mindist = dist;
	elseif dist > maxdist
		maxdist = dist;
	end
end

out.fes = counteval;
out.successrate = 1;
out.fmin = fmin;
out.stdx = std(x, [], 2);
out.k = K;
out.popsize = NP;
out.de_f = F;
out.fmean = mean(f);
out.cond = cond(x);
out.mindist = mindist;
out.maxdist = maxdist;

try
	while true
		
		% cluster
		c = cluster(linkage(x', 'ward', 'euclidean'), 'maxclust', round(K));
		c = de_sortcluster(f, round(K), c);
		
		if opts.plotting
			
			% draw contour and solutions
			if opts.de_more_plot
				subplot(4, 2, [1, 3]);
			else
				subplot(2, 2, [1, 3]);
			end
			
			hold off;
			[contourx, contoury, contourz] = ...
				drawcontour(fitfun, L, R, U, contourx, contoury, contourz);
			hold on;
			
			for i = 1 : length(c)
				selectcolor = mod(mod(i, length(c)), 6);
				if selectcolor == 0
					plot(x(1, c == i), x(2, c == i), 'kx', 'MarkerSize', 12, 'LineWidth', 3);
				elseif selectcolor == 1
					plot(x(1, c == i), x(2, c == i), 'bx', 'MarkerSize', 12, 'LineWidth', 3);
				elseif selectcolor == 2
					plot(x(1, c == i), x(2, c == i), 'yx', 'MarkerSize', 12, 'LineWidth', 3);
				elseif selectcolor == 3
					plot(x(1, c == i), x(2, c == i), 'cx', 'MarkerSize', 12, 'LineWidth', 3);
				elseif selectcolor == 4
					plot(x(1, c == i), x(2, c == i), 'mx', 'MarkerSize', 12, 'LineWidth', 3);
				else
					plot(x(1, c == i), x(2, c == i), 'gx', 'MarkerSize', 12, 'LineWidth', 3);
				end
			end
			
			% draw fitness
			if opts.de_more_plot
				subplot(4, 2, 2);
			else
				subplot(2, 2, 2);
			end
			
			hold off;
			semilogy(out.fes, out.fmin - ftarget, 'r');
			hold on;
			semilogy(out.fes, out.fmean - ftarget, 'b');
			xlabel('FEs');
			ylabel('fmin - ftarget');
			
			% draw standard deviation
			if opts.de_more_plot
				subplot(4, 2, 4);
			else
				subplot(2, 2, 4);
			end
			
			hold off;
			semilogy(out.fes, out.stdx(1, :));
			hold on;
			semilogy(out.fes, out.stdx(2, :));
			semilogy(out.fes, out.maxdist, 'r');
			semilogy(out.fes, out.mindist, 'r');
			xlabel('FEs');
			ylabel('STD');
			
			% more information
			if opts.de_more_plot
				subplot(4, 2, 5);
				plot(out.fes, out.k);
				xlabel('FEs');
				ylabel('K');
				subplot(4, 2, 6);
				plot(out.fes, out.successrate);
				xlabel('FEs');
				ylabel('SR');
				subplot(4, 2, 7);
				semilogy(out.fes, out.cond);
				xlabel('FEs');
				ylabel('Condition Number');
				subplot(4, 2, 8);
				semilogy(out.fes, out.de_f);
				xlabel('FEs');
				ylabel('F');
			end
			
			pause(0.02);
		end
		
		% termination criteria
		if counteval > opts.maxfunevals - NP
			stopflag = 'maxfunevals';
			break;
		elseif fmin <= ftarget
			stopflag = 'fitness';
			break;
		elseif std(f) <= 100 * eps * max(abs(f)) || std(f) <= 100 * realmin
			stopflag = 'TolFun';
			break;
		elseif countiter > opts.de_maxiter
			stopflag = 'maxiter';
			break;
		end
		
		% mutation/crossover/repair/selection
		countiter = countiter + 1;
		v = de_mutate(x, K, NP, F, c, opts);
		u = de_crossover(x, v, c, K, opts);
		u = de_repair(u, opts);
		[x, f, counteval, countsuccess] = de_select(x, u, f, counteval, fitfun, varargin{:});
		% 	[x, f, counteval] = de_regular(x, K, c, opts, f, counteval, fitfun, varargin{:});
		out = de_record(out, x, f, c, K, NP, F, counteval, countsuccess, opts);
		[f, idx] = sort(f);
		x = x(:, idx);
		xmin = x(:, 1);
		fmin = f(1);
		bestever.xmin = xmin;
		bestever.fmin = fmin;
		
		% adjust cluster number
		[~, prevNP] = size(x);
		[K, NP, F] = de_adjust(out, D, K, NP, F, counteval, opts);
		
		if NP < prevNP
			x = x(:, 1 : NP);
			f = f(1 : NP);
		end
	end
catch ME
	disp(getReport(ME));
	fprintf('NP = %f, K = %f, F = %f\n', NP, K, F);
	fprintf('St.D of f = %f', std(f));
	disp('Aborting the current evolution...');
end
end
