function f = cec15_f15(x)
% CEC15_F15 CEC'2015 function 15
% Composition Function 7 (N=10)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 15) - 1500;
end
