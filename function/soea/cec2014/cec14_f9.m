function f = cec14_f9(x)
% CEC14_F9 CEC'2014 function 9
% Shifted and Rotated Rastrigin¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 9) - 900;
end
