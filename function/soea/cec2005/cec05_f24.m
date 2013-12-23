function f = cec05_f24(x)
% CEC05_F24 CEC'2005 function 24
% Rotated Hybrid Composition Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 24) - 260;
if f < 1e-8
	f = 0;
end
end
