function f = cec05_f7(x)
% CEC05_F7 CEC'2005 function 7
% Shifted Rotated Griewank's Function without Bounds
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 7) + 180;
if f < 1e-8
	f = 0;
end
end
