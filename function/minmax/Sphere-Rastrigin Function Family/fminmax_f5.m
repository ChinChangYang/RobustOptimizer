function f = fminmax_f5(x, y)
% Sphere-rastrigin function (Type 1)
% Global optimum value: f(0.1, 0.1) = 0
%
% Property
% * xy-correlation: No
shift = 0.1;
if nargin == 0
	f = shift;
	return;
end

x = x - shift;
y = y - shift;
X = x;
Y = y;
a = 1;
b = -2;
f = a * sphere(X) + b * rastrigin(Y);
end
