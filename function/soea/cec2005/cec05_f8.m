function f = cec05_f8(x)
% CEC05_F8 CEC'2005 function 8
% Shifted Rotated Ackley's Function with Global Optimum on Bounds
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 8) + 140;
if f < 1e-8
	f = 0;
end
end
