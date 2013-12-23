function f = cloppower(x)
%CLOPPOWER A test function with POWER definition in CLOP paper
%The function values are Bernoulli distributed with a parameter p defined
%as f(x) = 1 / (1 + exp(-r(x))), where r(x) is defined as
%0.05(x+1)^2-((x+1)/2)^20.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
r = sum(0.05 * (x + 1).^2 - ((x + 1) / 2).^20);
p = 1 / (1 + exp(-r));

if rand < p
	f = 0;
else
	f = 1;
end
end
