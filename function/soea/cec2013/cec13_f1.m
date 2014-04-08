function f = cec13_f1(x)
% CEC13_F1 CEC'2013 function 1
% Sphere function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 1) + 1400;
end
