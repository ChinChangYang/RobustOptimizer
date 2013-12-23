function f = maxminmax_f28(x,y,z)
% Sphere-sphere-rastrigin function (Type 4)
% Global optimum value: f(0.1,0.1,0.1) = 0
%
% Property
% * Multimodal
% * xy-correlation: No
% * xz-correlation: No
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
Y = y;
Z = y + z;
a = 1;
b = 1;
c = 2;
f = - a * sum(X.^2) ...
	+ b * sum(Y.^2) ...
	- c * rastrigin(Z);
end
