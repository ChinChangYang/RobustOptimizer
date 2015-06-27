function f = cec15_f11(x)
% CEC15_F11 CEC'2015 function 11
% Composition Function 3 (N=5)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 11) - 1100;
end
