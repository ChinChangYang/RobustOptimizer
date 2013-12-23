function f = cec05_f2(x)
% CEC05_F2 CEC'2005 function 2
% Shifted Schwefel's Problem 1.2
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 2) + 450;
if f < 1e-8
	f = 0;
end
end
