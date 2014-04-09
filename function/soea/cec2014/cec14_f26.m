function f = cec14_f26(x)
% CEC14_F26 CEC'2014 function 26
% Composition Function 4
% N = 5
% sigma = [10, 10, 10, 10, 10]
% lambda = [ 0.25, 1, 1e-7, 2.5, 10]
% bias = [0, 100, 200, 300, 400]
% g1: Rotated Schwefel's Function F11・
% g2: Rotated HappyCat Function F13・
% g3: Rotated High Conditioned Elliptic Function F1・
% g4: Rotated Weierstrass Function F6・
% g5: Rotated Griewank・s Function F7・
x = reshape(x, numel(x), 1);
f = cec14_func(x, 26) - 2600;
end
