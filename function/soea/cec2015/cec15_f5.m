function f = cec15_f5(x)
% CEC15_F5 CEC'2015 function 5
% Shifted and Rotated Schwefel¡¦s Function
x = reshape(x, numel(x), 1);
f = cec15_func(x, 5) - 500;
end
