function f = cec14_f3(x)
% CEC14_F3 CEC'2014 function 3
% Rotated Discus Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 3) - 300;
end
