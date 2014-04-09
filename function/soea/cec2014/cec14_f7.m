function f = cec14_f7(x)
% CEC14_F7 CEC'2014 function 7
% Shifted and Rotated Griewank¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 7) - 700;
end
