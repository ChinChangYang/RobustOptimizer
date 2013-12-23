function f = cec05_f10(x)
% CEC05_F10 CEC'2005 function 10
% Shifted Rotated Rastrigin's Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 10) + 330;
if f < 1e-8
	f = 0;
end
end
