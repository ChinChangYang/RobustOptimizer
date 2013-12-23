function f = cec05_f22(x)
% CEC05_F22 CEC'2005 function 22
% Rotated Hybrid Composition Function with High Condition Number Matrix
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 22) - 360;
if f < 1e-8
	f = 0;
end
end
