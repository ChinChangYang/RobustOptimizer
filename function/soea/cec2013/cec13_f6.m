function f = cec13_f6(x)
% CEC13_F6 CEC'2013 function 6
% Rotated Rosenbrock's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 6) + 900;
end

