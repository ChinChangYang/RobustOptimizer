function f = cec14_f14(x)
% CEC14_F14 CEC'2014 function 14
% Shifted and Rotated HGBat Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 14) - 1400;
end
