function f = cec15_f3(x)
% CEC15_F3 CEC'2015 function 3
% Shifted and Rotated Ackley¡¦s Function
x = reshape(x, numel(x), 1);
f = cec15_func(x, 3) - 300;
end
