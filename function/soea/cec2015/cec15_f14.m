function f = cec15_f14(x)
% CEC15_F14 CEC'2015 function 14
% Composition Function 6 (N=7)
x = reshape(x, numel(x), 1);
f = cec15_func(x, 14) - 1400;
end
