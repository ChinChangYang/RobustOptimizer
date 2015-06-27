function f = cec15_f6(x)
% CEC15_F6 CEC'2015 function 6
% Hybrid Function 1 (N=3)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 6) - 600;
end
