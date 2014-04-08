function f = cec13_f19(x)
% CEC13_F19 CEC'2013 function 19
% Expanded Griewank's plus Rosenbrock's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 19) - 500;
end

