function f = cec14_f28(x)
% CEC14_F28 CEC'2014 function 28
% Composition Function 6
% N = 5
% sigma = [10, 20, 30, 40, 50]
% lambda = [ 2.5, 10, 2.5, 5e-4,1e-6]
% bias = [0, 100, 200, 300, 400]
% g1: Rotated Expanded Griewank・s plus Rosenbrock・s Function F15・
% g2: Rotated HappyCat Function F13・
% g3: Rotated Schwefel's Function F11・
% g4: Rotated Expanded Scaffer・s F6 Function F16・
% g5: Rotated High Conditioned Elliptic Function F1・
x = reshape(x, numel(x), 1);
f = cec14_func(x, 28) - 2800;
end
