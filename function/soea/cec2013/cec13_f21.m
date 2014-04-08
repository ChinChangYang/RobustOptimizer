function f = cec13_f21(x)
% CEC13_F21 CEC'2013 function 21
% Composition Function 1
% n = 5
% g1: Rotated Rosenbrock's Function f6
% g2: Rotated High Conditioned Elliptic Function f2
% g3: Rotated Bent Cigar Function f4
% g4: Rotated Discus Function f3
% g5: Sphere Function f1
x = reshape(x, numel(x), 1);
f = cec13_func(x, 21) - 700;
end

