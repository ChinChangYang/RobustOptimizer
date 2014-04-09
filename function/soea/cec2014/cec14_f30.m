function f = cec14_f30(x)
% CEC14_F30 CEC'2014 function 30
% Composition Function 8
% N = 3
% sigma = [10, 30, 50]
% lambda = [1, 1, 1]
% bias = [0, 100, 200]
% g1: Hybrid Function 4 F20¡¦
% g2: Hybrid Function 5 F21¡¦
% g3: Hybrid Function 6 F22¡¦
x = reshape(x, numel(x), 1);
f = cec14_func(x, 30) - 3000;
end
