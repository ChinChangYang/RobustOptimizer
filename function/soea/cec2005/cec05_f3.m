function f = cec05_f3(x)
% CEC05_F3 CEC'2005 function 3
% Shifted Rotated High Conditioned Elliptic Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 3) + 450;
if f < 1e-8
	f = 0;
end
end
