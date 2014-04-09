function f = cec14_f23(x)
% CEC14_F23 CEC'2014 function 23
% Composition Function 1
% N= 5, sigma = [10, 20, 30, 40, 50]
% lambda = [ 1, 1e-6, 1e-26, 1e-6, 1e-6]
% bias = [0, 100, 200, 300, 400]
% g1: Rotated Rosenbrock・s Function F4・
% g2: High Conditioned Elliptic Function F1・
% g3 Rotated Bent Cigar Function F2・
% g4: Rotated Discus Function F3・
% g5: High Conditioned Elliptic Function F1・
x = reshape(x, numel(x), 1);
f = cec14_func(x, 23) - 2300;
end
