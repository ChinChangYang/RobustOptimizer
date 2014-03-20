function f = cec14_f5(x)
% CEC14_F5 CEC'2014 function 5
% Shifted and Rotated Ackley¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 5) - 500;
if f < 1e-8
	f = 0;
end
end
