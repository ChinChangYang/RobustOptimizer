function [xmin, ...      % minimum search point of last iteration
	fmin, ...      % function value of xmin
	counteval, ... % number of function evaluations done
	stopflag, ...  % stop criterion reached
	out, ...     % struct with various histories and solutions
	bestever ... % struct containing overall best solution (for convenience)
	] = desacs( ...
	fitfun, ...    % name of objective/fitness function
	xstart, ...    % objective variables initial point, determines N
	insteps, ...   % initial coordinate wise population steps
	inopts, ...    % options struct, see defopts below
	varargin )     % arguments passed to objective function
%DESaCS Differential evolution with self-adaptive cluster size

%% 檢查輸入的參數
if nargin < 1
	%% 自我測試
	fitfun = 'bbob12_f15';
	xstart = zeros(2, 1);
	insteps = 10 * ones(2, 1);
	inopts.maxfunevals = 1e6;
	inopts.plotting = true;
	inopts.de_k = 3;
	inopts.de_globalpopsize = 150;
	varargin = [];
	
	s = RandStream('mt19937ar', 'Seed', 14);
	RandStream.setGlobalStream(s);
	
	inopts.ftarget = feval(fitfun, 'xopt') + 1e-8;
	[xmin, fmin, counteval, stopflag, out, bestever] = ...
		desacs(fitfun, xstart, insteps, inopts, varargin);
	
	if bestever.fmin <= inopts.ftarget
		disp('Success!');
	else
		disp('Failed!');
	end
	
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

%% 初始化輸出資料
out.kmeansde.counteval = 0;
out.kmeansde.k = 0;
out.kmeansde.popsize = 0;
out.hugede.counteval = 0;
out.hugede.k = 0;
out.hugede.popsize = 0;

%% 執行小人口數的 DESaCS
kmeansdeopts = opts;
kmeansdeopts.maxfunevals = opts.maxfunevals;

[xmin, fmin, kmeansdecounteval, stopflag, kmeansdeout, kmeansdebestever] = ...
	kmeansde(fitfun, xstart, insteps, kmeansdeopts, varargin{:});

counteval = kmeansdecounteval;
out.fes = kmeansdeout.fes;
out.fmin = kmeansdeout.fmin;
out.kmeansde.counteval = kmeansdecounteval;
out.kmeansde.k = kmeansdeout.k(end);
out.kmeansde.popsize = kmeansdeopts.de_globalpopsize;

bestever.xmin = kmeansdebestever.xmin;
bestever.fmin = kmeansdebestever.fmin;

if strcmp('fitness', stopflag)
	stopflag = 'normal desacs fitness';
	return;
end

if strcmp('maxfunevals', stopflag)
	stopflag = 'normal desacs maxfunevals';
	return;
end

%% 執行大人口數的 DESaCS
kmeansdeopts = opts;
kmeansdeopts.maxfunevals = opts.maxfunevals - counteval;
kmeansdeopts.de_globalpopsize = min(5000, max(3.0, opts.de_globalpopsize * (opts.maxfunevals - counteval) / out.kmeansde.counteval));
kmeansdeopts.de_maxiter = Inf;

[xmin, fmin, kmeansdecounteval, stopflag, kmeansdeout, kmeansdebestever] = ...
	kmeansde(fitfun, xstart, insteps, kmeansdeopts, varargin{:});

counteval = counteval + kmeansdecounteval;
kmeansdeout.fes = kmeansdeout.fes + out.fes(end);
out.fes = [out.fes, kmeansdeout.fes];
kmeansdeout.fmin(kmeansdeout.fmin > out.fmin(end)) = out.fmin(end);
out.fmin = [out.fmin, kmeansdeout.fmin];
out.hugede.counteval = kmeansdecounteval;
out.hugede.k = kmeansdeout.k(end);
out.hugede.popsize = kmeansdeopts.de_globalpopsize;

if kmeansdebestever.fmin <= bestever.fmin
	bestever.xmin = kmeansdebestever.xmin;
	bestever.fmin = kmeansdebestever.fmin;
end

if strcmp('fitness', stopflag)
	stopflag = 'huge desacs fitness';
	return;
end

if strcmp('maxfunevals', stopflag)
	stopflag = 'huge desacs maxfunevals';
	return;
end

if strcmp('TolFun', stopflag)
	stopflag = 'huge desacs TolFun';
	return
end

if strcmp('maxiter', stopflag)
	stopflag = 'huge desacs maxiter';
	return;
end

