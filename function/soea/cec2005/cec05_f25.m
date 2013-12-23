function f = cec05_f25(x)
% CEC05_F25 CEC'2005 function 25
% Rotated Hybrid Composition Function without bounds
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 25) - 260;
if f < 1e-8
	f = 0;
end
end
