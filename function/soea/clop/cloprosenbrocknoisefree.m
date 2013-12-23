function f = cloprosenbrocknoisefree(x)
%CLOPROSENBROCKNOISEFREE A test function with ROSENBROCK definition in CLOP paper
%The function values are defined as f(x) = 1 / (1 + exp(-r(x))), where r(x)
%is defined as 1.0-0.1((1-a)^2+(b-a^2)^2), where a = 4x_1 and b = 10x_2+4.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
%
% Minimizer: [0.25; -0.3]
% Minimum function value: 2.689414213699951e-01
a = 4 * x(1);
b = 10 * x(2) + 4;
r = 1.0 - 0.1 * ((1 - a)^2 + (b - a^2)^2);
f = 1 - 1 / (1 + exp(-r));
end
