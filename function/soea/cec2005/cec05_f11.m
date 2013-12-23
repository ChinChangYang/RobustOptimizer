function f = cec05_f11(x)
% CEC05_F11 CEC'2005 function 11
% Shifted Rotated Weierstrass Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 11) - 90;
if f < 1e-8
	f = 0;
end
end
