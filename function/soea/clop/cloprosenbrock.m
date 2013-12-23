function f = cloprosenbrock(x)
%CLOPROSENBROCK A test function with ROSENBROCK definition in CLOP paper
%The function values are Bernoulli distributed with a parameter p defined
%as f(x) = 1 / (1 + exp(-r(x))), where r(x) is defined as
%1.0-0.1((1-a)^2+(b-a^2)^2), where a = 4x_1 and b = 10x_2+4.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
a = 4 * x(1);
b = 10 * x(2) + 4;
r = 1.0 - 0.1 * ((1 - a)^2 + (b - a^2)^2);
p = 1 / (1 + exp(-r));

if rand < p
	f = 0;
else
	f = 1;
end
end
