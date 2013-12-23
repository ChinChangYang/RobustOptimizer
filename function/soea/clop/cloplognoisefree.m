function f = cloplognoisefree(x)
%CLOPLOGNOISEFREE A test function with LOG definition in CLOP paper
%The function values are defined as f(x) = 1 / (1 + exp(-r(x))), where r(x)
%is defined as 2log(4x+4.1)-4x-3.
%
% Note that the CLOP paper defines the test function for maximization. In
% this project, clopangle defines an invert version of the test function
% for minimization.
%
% Minimizer: -0.52500
% Minimum 1D function value: 3.807669091787934e-01
% Minimum 5D function value: 8.080408170692532e-02
r = sum(2 * log(4 * x + 4.1) - 4 * x - 3);
f = 1 - 1 / (1 + exp(-r));
end
