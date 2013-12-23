function f = dampedcosine(x,y)
% Damped cosine wave function
% Global min-max value #1: f(40.8829, 100) = ?
% Global min-max value #2: f(40.8829, -100) = ?
if nargin == 0
	f = [40.8829,	40.8829; ...
		100,		-100];
	
	return;
end

x = 0.05 * x + 5;
y = 0.05 * y + 5;
f = cos(sqrt(sum(x.^2) + sum(y.^2))) / (sqrt(sum(x.^2) + sum(y.^2)) + 10);
end
