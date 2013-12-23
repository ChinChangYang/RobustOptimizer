function f = cec05_f19(x)
% CEC05_F19 CEC'2005 function 19
% Rotated Hybrid Composition Function with narrow basin global optimum
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 19) - 10;
if f < 1e-8
	f = 0;
end
end
