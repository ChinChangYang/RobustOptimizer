function f = cec15_f10(x)
% CEC15_F10 CEC'2015 function 10
% Composition Function 2 (N=3)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 10) - 1000;
end
