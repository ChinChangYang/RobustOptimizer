function f = cec05_f15(x)
% CEC05_F15 CEC'2005 function 15
% Shifted sphere function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 15) - 120;
if f < 1e-8
	f = 0;
end
end
