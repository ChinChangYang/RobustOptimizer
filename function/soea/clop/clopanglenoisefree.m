function f = clopanglenoisefree(x)
%CLOPANGLENOISEFREE A test function with ANGLE definition in CLOP paper
%The function values are Bernoulli distributed with a parameter p defined
%as f(x) = 1 / (1 + exp(-r(x))), where r(x) is defined as
%the ANGLE function.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
%
% Minimizer: -0.2
% Minimum function value: 2.689414213699976e-01
r = 0;
for i = 1 : numel(x)
	if x(i) < -0.2
		r = r + 1 + sqrt(2) - 2 * sqrt(0.3 - x(i));
	else
		r = r + 1 + sqrt(2) - sqrt(x(i) + 2.2);
	end
end

f = 1 - 1 / (1 + exp(-r));
end
