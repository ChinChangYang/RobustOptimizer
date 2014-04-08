function f = cec13_f16(x)
% CEC13_F16 CEC'2013 function 16
% Rotated Katsuura Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 16) - 200;
end

