function f = cec13_f10(x)
% CEC13_F10 CEC'2013 function 10
% Rotated Griewank's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 10) + 500;
end

