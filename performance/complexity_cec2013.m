function [ret, t0, t1, t2] = complexity_cec2013(solver, D)
%COMPLEXITY_CEC2013 Algorithm Complexity defined in CEC'2013

% Measure T0
tic;
for i = 1 : 1000000
	x = 5.55;
	x = x + x; 
	x = x ./ 2; 
	x = x * x; 
	x = sqrt(x); 
	x = log(x); 
	x = exp(x); 
	y = x / x; %#ok<NASGU>
end
t0 = toc;

% Measure T1
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);
maxfunevals = 200000;
tic;
feval(solver, 'cec13_f14', lb, ub, maxfunevals);
t1 = toc;

% Measure T2
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);
maxfunevals = 200000;
all_t2 = zeros(1, 5);
for i = 1 : 5
	tic;
	feval(solver, 'cec13_f3', lb, ub, maxfunevals);
	all_t2(i) = toc;
end
t2 = mean(all_t2);

% Measure complexity of the solver
ret = (t2 - t1) / t0;
end

