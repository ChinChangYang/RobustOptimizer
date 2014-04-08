function f = cec13_f28(x)
% CEC13_F28 CEC'2013 function 28
% Composition Function 8
% n = 5
% g1: Rotated Expanded Griewank's plus Rosenbrock's Function f19
% g2: Rotated Schaffers F7 Function f7
% g3: Rotated Schwefel's Function f15
% g4: Rotated Expanded Scaffer's F6 Function f20
% g5: Sphere Function f1
x = reshape(x, numel(x), 1);
f = cec13_func(x, 28) - 1400;
end

