function f = cec15_f1(x)
% CEC15_F1 CEC'2015 function 1
% Rotated High Conditioned Elliptic Function
x = reshape(x, numel(x), 1);
f = cec15_func(x, 1) - 100;
end
