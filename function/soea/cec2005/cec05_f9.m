function f = cec05_f9(x)
% CEC05_F9 CEC'2005 function 9
% Shifted Rastrigin's Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 9) + 330;
if f < 1e-8
	f = 0;
end
end
