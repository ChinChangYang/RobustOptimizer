function f = ellipsoidalD2d45(x)
% ELLIPSOIDALD2d45 2-D Ellipsoidal function rotated by 45 degree
%
% * unimodal
% * conditioning is 1e6
%
M = numel(x);
x = reshape(x, 1, M);
B = [cos(pi/4), -sin(pi/4); sin(pi/4), cos(pi/4)];
x = x * B(1:M, 1:M);
f = sum(10.^(6/(M-1)*(0:M-1)) .* x .^2);
end

