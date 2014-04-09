function f = cec14_f25(x)
% CEC14_F25 CEC'2014 function 25
% Composition Function 3
% N = 3
% sigma = [10, 30, 50]
% lambda = [0.25, 1, 1e-7]
% bias = [0, 100, 200]
% g1: Rotated Schwefel's Function F11・
% g2: Rotated Rastrigin・s Function F9・
% g3: Rotated High Conditioned Elliptic Function F1・
x = reshape(x, numel(x), 1);
f = cec14_func(x, 25) - 2500;
end
