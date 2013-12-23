function f = cec13_f9(x)
% CEC13_F9 CEC'2013 function 9
% Rotated Weierstrass Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 9) + 600;
if f < 1e-8
	f = 0;
end
end

