function f = clopstep(x)
%CLOPSTEP A test function with STEP definition in CLOP paper
%The function values are Bernoulli distributed with a parameter p defined
%as f(x) = 1 / (1 + exp(-r(x))), where r(x) is defined as
%the STEP function.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
r = 0;
for i = 1 : numel(x)
	if x(i) < -0.8
		r = r - 2;
	elseif x(i) < -0.3
		r = r - 2 + 6 * (x(i) + 0.8);
	elseif x(i) < 0.8
		r = r - (x(i) + 0.3) / 1.1;
	else
		r = r -2;
	end
end

p = 1 / (1 + exp(-r));

if rand < p
	f = 0;
else
	f = 1;
end
end
