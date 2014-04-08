function f = cec13_f11(x)
% CEC13_F11 CEC'2013 function 11
% Rastrigin's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 11) + 400;
end

