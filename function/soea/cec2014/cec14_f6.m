function f = cec14_f6(x)
% CEC14_F6 CEC'2014 function 6
% Shifted and Rotated Weierstrass Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 6) - 600;
if f < 1e-8
	f = 0;
end
end
