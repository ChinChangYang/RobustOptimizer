function f = clopflat(x)
%CLOPFLAT A test function with FLAT definition in CLOP paper
%The function values are Bernoulli distributed with a parameter p defined
%as f(x) = 1 / (1 + exp(-r(x))), where r(x) is defined as
%0.2/(1+6(x+0.6)^2+(x+0.6)^3).
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
r = sum(0.2 ./ (1 + 6 * (x + 0.6).^2 + (x + 0.6).^3));
p = 1 / (1 + exp(-r));

if rand < p
	f = 0;
else
	f = 1;
end
end
