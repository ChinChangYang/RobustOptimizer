function f = cec14_f15(x)
% CEC14_F15 CEC'2014 function 15
% Shifted and Rotated Expanded Griewank¡¦s plus Rosenbrock¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 15) - 1500;
end
