function f = maxminmax_f2(x,y,z)
% Sphere-sphere-sphere function (Type 2)
% Global max-min-max value: f(0.1,0.1,0.1) = 0
%
% Property
% * Unimodal
% * xy-correlation: Yes
% * xz-correlation: No
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
Z = z;
f = -sum(X.^2) ...
	+ sum(Y.^2) ...
	- sum(Z.^2);
end
