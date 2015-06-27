function f = cec15_f8(x)
% CEC15_F8 CEC'2015 function 8
% Hybrid Function 3(N=5)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 8) - 800;
end
