function f = cec13_f15(x)
% CEC13_F15 CEC'2013 function 15
% Rotated Schwefel's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 15) - 100;
if f < 1e-8
	f = 0;
end
end

