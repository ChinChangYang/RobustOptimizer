function f = cec14_f2(x)
% CEC14_F2 CEC'2014 function 2
% Rotated Bent Cigar Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 2) - 200;
end
