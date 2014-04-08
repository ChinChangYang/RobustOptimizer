function f = cec13_f13(x)
% CEC13_F13 CEC'2013 function 13
% Non-Continuous Rotated Rastrigin's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 13) + 200;
end

