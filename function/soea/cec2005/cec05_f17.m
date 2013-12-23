function f = cec05_f17(x)
% CEC05_F17 CEC'2005 function 17
% Shifted sphere function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 17) - 120;
if f < 1e-8
	f = 0;
end
end
