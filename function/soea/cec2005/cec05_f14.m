function f = cec05_f14(x)
% CEC05_F14 CEC'2005 function 14
% Shifted Rotated Expanded Scaffer's F6 Function
x = reshape(x, numel(x), 1);
f = benchmark_func(x, 14) + 300;
if f < 1e-8
	f = 0;
end
end
