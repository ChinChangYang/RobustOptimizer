function f = cec13_f4(x)
% CEC13_F4 CEC'2013 function 4
% Rotated Discus Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 4) + 1100;
if f < 1e-8
	f = 0;
end
end

