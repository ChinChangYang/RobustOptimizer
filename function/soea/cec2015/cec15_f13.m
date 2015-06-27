function f = cec15_f13(x)
% CEC15_F13 CEC'2015 function 13
% Composition Function 5 (N=5)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 13) - 1300;
end
