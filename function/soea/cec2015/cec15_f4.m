function f = cec15_f4(x)
% CEC15_F4 CEC'2015 function 4
% Shifted and Rotated Rastrigin¡¦s Function
x = reshape(x, numel(x), 1);
f = cec15_func(x, 4) - 400;
end
