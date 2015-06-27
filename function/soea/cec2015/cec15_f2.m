function f = cec15_f2(x)
% CEC15_F2 CEC'2015 function 2
% Rotated Cigar Function
x = reshape(x, numel(x), 1);
f = cec15_func(x, 2) - 200;
end
