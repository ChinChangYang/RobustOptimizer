function f = cloppowernoisefree(x)
%CLOPPOWERNOISEFREE A test function with POWER definition in CLOP paper
%The function values are Bernoulli distributed with a parameter p defined
%as f(x) = 1 / (1 + exp(-r(x))), where r(x) is defined as
%0.05(x+1)^2-((x+1)/2)^20.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
%
% Minimizer: 6.093213380198702e-01
% Minimum function value: 4.708963893093869e-01
r = sum(0.05 * (x + 1).^2 - ((x + 1) / 2).^20);
f = 1 - 1 / (1 + exp(-r));
end
