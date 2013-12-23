function f = maxminmax_f45(x,y,z)
% Rastrigin-sphere-rastrigin function (Type 5)
% Global optimum value: f(0.1,0.1,0.1) = 0
%
% Property
% * Multimodal
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
f = - a * rastrigin(X) ...
	+ b * sum(Y.^2) ...
	- c * rastrigin(Z);
end
