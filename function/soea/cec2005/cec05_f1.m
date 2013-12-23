function f = cec05_f1(x)
% CEC05_F1 CEC'2005 function 1
% Shifted sphere function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 1) + 450;
if f < 1e-8
	f = 0;
end
end
