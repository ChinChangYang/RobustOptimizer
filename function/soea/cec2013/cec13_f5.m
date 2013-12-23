function f = cec13_f5(x)
% CEC13_F5 CEC'2013 function 5
% Different Powers Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 5) + 1000;
if f < 1e-8
	f = 0;
end
end

