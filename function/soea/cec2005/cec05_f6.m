function f = cec05_f6(x)
% CEC05_F6 CEC'2005 function 6
% Shifted Rosenbrock's Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 6) - 390;
if f < 1e-8
	f = 0;
end
end
