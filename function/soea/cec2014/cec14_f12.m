function f = cec14_f12(x)
% CEC14_F12 CEC'2014 function 12
% Shifted and Rotated Katsuura Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 12) - 1200;
if f < 1e-8
	f = 0;
end
end
