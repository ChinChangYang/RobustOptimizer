function f = generic_solver_err(x, varargin)
%GENERIC_SOLVER_ERR Generic solver error
% x is a solver parameter vector
solver	= varargin{1};
D		= varargin{2};
fitfun	= varargin{3};
maxfes	= varargin{4};
lb		= -100 * ones(D, 1);
ub		= 100 * ones(D, 1);

[~, T] = size(x);
f = zeros(1, T);

for t = 1 : T
	options = setDEoptions(x, t);
	
	[~, f(t), ~] = ...
		feval(solver, fitfun, lb, ub, maxfes, options);
end
end

