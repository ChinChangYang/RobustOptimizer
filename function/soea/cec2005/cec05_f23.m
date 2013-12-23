function f = cec05_f23(x)
% CEC05_F23 CEC'2005 function 23
% Non-Continuous Rotated Hybrid Composition Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 23) - 360;
if f < 1e-8
	f = 0;
end
end
