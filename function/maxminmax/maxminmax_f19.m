function f = maxminmax_f19(x,y,z)
% Sphere-rastrigin-sphere function (Type 3)
% Global max-min-max value: f(0.1,0.1,0.1) = 0
%
% Property
% * Multimodal
% * xy-correlation: No
% * xz-correlation: Yes
% * yz-correlation: No
shift = 0.1;
if nargin == 0
	f = shift;
	return;
end

x = x - shift;
y = y - shift;
z = z - shift;
X = x;
Y = y;
Z = x + z;
a = 1;
b = 1;
c = 1;
f = - a * sum(X.^2) ...
	+ b * rastrigin(Y) ...
	- c * sum(Z.^2);
end
