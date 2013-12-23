function f = cec05_f20(x)
% CEC05_F20 CEC'2005 function 20
% Rotated Hybrid Composition Function with Global Optimum on the Bounds
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 20) - 10;
if f < 1e-8
	f = 0;
end
end
