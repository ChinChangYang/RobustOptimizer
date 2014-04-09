function f = cec14_f27(x)
% CEC14_F27 CEC'2014 function 27
% Composition Function 5
% N = 5
% sigma = [10, 10, 10, 20, 20]
% lambda = [10, 10, 2.5, 25, 1e-6]
% bias = [0, 100, 200, 300, 400]
% g1: Rotated HGBat Function F14・
% g2: Rotated Rastrigin・s Function F9・
% g3: Rotated Schwefel's Function F11・
% g4: Rotated Weierstrass Function F6・
% g5: Rotated High Conditioned Elliptic Function F1・
x = reshape(x, numel(x), 1);
f = cec14_func(x, 27) - 2700;
end
