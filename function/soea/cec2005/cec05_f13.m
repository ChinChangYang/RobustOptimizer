function f = cec05_f13(x)
% CEC05_F13 CEC'2005 function 13
% Shifted Expanded Griewank's plus Rosenbrock's Function (F8F2)
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 13) + 130;
if f < 1e-8
	f = 0;
end
end
