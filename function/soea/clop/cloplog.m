function f = cloplog(x)
%CLOPLOG A test function with LOG definition in CLOP paper
%The function values are Bernoulli distributed with a parameter p defined
%as f(x) = 1 / (1 + exp(-r(x))), where r(x) is defined as
%2log(4x+4.1)-4x-3.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
r = sum(2 * log(4 * x + 4.1) - 4 * x - 3);
p = 1 / (1 + exp(-r));

if rand < p
	f = 0;
else
	f = 1;
end
end
