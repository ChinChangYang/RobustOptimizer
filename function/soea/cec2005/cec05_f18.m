function f = cec05_f18(x)
% CEC05_F18 CEC'2005 function 18
% Rotated Hybrid Composition Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 18) - 10;
if f < 1e-8
	f = 0;
end
end
