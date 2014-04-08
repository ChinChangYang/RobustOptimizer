function f = cec13_f12(x)
% CEC13_F12 CEC'2013 function 12
% Rotated Rastrigin's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 12) + 300;
end

