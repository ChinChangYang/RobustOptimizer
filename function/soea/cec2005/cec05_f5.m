function f = cec05_f5(x)
% CEC05_F5 CEC'2005 function 5
% Schwefel's Problem 2.6 with Global Optimum on Bounds
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 5) + 310;
if f < 1e-8
	f = 0;
end
end
