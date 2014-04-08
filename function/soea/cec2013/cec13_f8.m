function f = cec13_f8(x)
% CEC13_F8 CEC'2013 function 8
% Rotated Ackley's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 8) + 700;
end

