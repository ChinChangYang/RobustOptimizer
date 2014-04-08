function f = cec13_f24(x)
% CEC13_F24 CEC'2013 function 24
% Composition Function 4
% n = 3
% g1: Rotated Schwefel's Function f15
% g2: Rotated Rastrigin's Function f12
% g3: Rotated Weierstrass Function f9
x = reshape(x, numel(x), 1);
f = cec13_func(x, 24) - 1000;
end

