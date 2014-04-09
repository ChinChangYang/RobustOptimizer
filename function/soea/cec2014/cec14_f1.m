function f = cec14_f1(x)
% CEC14_F1 CEC'2014 function 1
% Rotated High Conditioned Elliptic Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 1) - 100;
end
