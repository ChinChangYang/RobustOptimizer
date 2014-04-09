function f = cec14_f29(x)
% CEC14_F29 CEC'2014 function 29
% Composition Function 7
% N = 3
% sigma = [10, 30, 50]
% lambda = [1, 1, 1]
% bias = [0, 100, 200]
% g1: Hybrid Function 1 F17¡¦
% g2: Hybrid Function 2 F18¡¦
% g3: Hybrid Function 3 F19¡¦
x = reshape(x, numel(x), 1);
f = cec14_func(x, 29) - 2900;
end
