function f = cec13_f18(x)
% CEC13_F18 CEC'2013 function 18
% Rotated Lunacek Bi_rastrigin Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 18) - 400;
if f < 1e-8
	f = 0;
end
end

