function f = cec15_f7(x)
% CEC15_F7 CEC'2015 function 7
% Hybrid Function 2 (N=4)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 7) - 700;
end
