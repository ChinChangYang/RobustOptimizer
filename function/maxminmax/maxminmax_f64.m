function f = maxminmax_f64(x,y,z)
% Rastrigin-rastrigin-rastrigin function (Type 8)
% Global optimum value: f(0.1,0.1,0.1) = 0
%
% Property
% * Multimodal
% * xy-correlation: Yes
% * xz-correlation: Yes
% * yz-correlation: Yes
shift = 0.1;
if nargin == 0
	f = shift;
	return;
end

x = x - shift;
y = y - shift;
z = z - shift;
X = x;
Y = x + y;
Z = y + z;
a = 1;
b = 2;
c = 2;
f = - a * rastrigin(X) ...
	+ b * rastrigin(Y) ...
	- c * rastrigin(Z);
end
