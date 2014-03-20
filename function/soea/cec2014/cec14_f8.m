function f = cec14_f8(x)
% CEC14_F8 CEC'2014 function 8
% Shifted Rastrigin¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 8) - 800;
if f < 1e-8
	f = 0;
end
end
