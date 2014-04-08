function f = cec13_f17(x)
% CEC13_F17 CEC'2013 function 17
% Lunacek Bi_rastrigin Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 17) - 300;
end

