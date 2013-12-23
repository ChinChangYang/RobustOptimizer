function f = cec13_f2(x)
% CEC13_F2 CEC'2013 function 2
% Rotated High Conditioned Elliptic Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 2) + 1300;
if f < 1e-8
	f = 0;
end
end

