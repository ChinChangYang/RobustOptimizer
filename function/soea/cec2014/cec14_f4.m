function f = cec14_f4(x)
% CEC14_F4 CEC'2014 function 4
% Shifted and Rotated Rosenbrock¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 4) - 400;
end
