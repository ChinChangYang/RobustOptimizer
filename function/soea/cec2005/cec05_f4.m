function f = cec05_f4(x)
% CEC05_F4 CEC'2005 function 4
% Shifted Schwefel's Problem 1.2 with Noise in Fitness
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 4) + 450;
if f < 1e-8
	f = 0;
end
end
