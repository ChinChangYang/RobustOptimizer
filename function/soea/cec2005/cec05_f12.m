function f = cec05_f12(x)
% CEC05_F12 CEC'2005 function 12
% Schwefel's Problem 2.13
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 12) + 460;
if f < 1e-8
	f = 0;
end
end
