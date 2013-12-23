function f = maxminmax_f5(x,y,z)
% Sphere-sphere-sphere function (Type 5)
% Global max-min-max value: f(0.1,0.1,0.1) = 0
%
% Property
% * Unimodal
% * xy-correlation: Yes
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
Y = x + y;
Z = x + z;
a = 1;
b = 1;
c = 2;
f = - a * sum(X.^2) ...
	+ b * sum(Y.^2) ...
	- c * sum(Z.^2);
end
