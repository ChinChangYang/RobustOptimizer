function f = cec05_f16(x)
% CEC05_F16 CEC'2005 function 16
% Shifted sphere function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 16) - 120;
if f < 1e-8
	f = 0;
end
end
