function f = cec15_f9(x)
% CEC15_F9 CEC'2015 function 9
% Composition Function 1 (N=3)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 9) - 900;
end
