function f = dampedsinus(x, y)
% Damped sinus function
% Global min-max value: f(100,-57.4863) = ?
x = 0.05 * x + 5;
y = 0.05 * y + 5;
f = sum(sin(x-y)) / sum(sqrt(x.^2 + y.^2));
end
