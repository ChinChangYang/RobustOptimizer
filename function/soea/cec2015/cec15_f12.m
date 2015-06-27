function f = cec15_f12(x)
% CEC15_F12 CEC'2015 function 12
% Composition Function 4 (N=5)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 12) - 1200;
end
