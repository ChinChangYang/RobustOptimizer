function f = cec14_f11(x)
% CEC14_F11 CEC'2014 function 11
% Shifted and Rotated Schwefel¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 11) - 1100;
if f < 1e-8
	f = 0;
end
end
