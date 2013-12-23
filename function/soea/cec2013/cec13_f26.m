function f = cec13_f26(x)
% CEC13_F26 CEC'2013 function 26
% Composition Function 6
% n = 5
% g1: Rotated Schwefel's Function f15
% g2: Rotated Rastrigin's Function f12
% g3: Rotated High Conditioned Elliptic Function f2
% g4: Rotated Weierstrass Function f9
% g5: Rotated Griewank's Function f10
x = reshape(x, numel(x), 1);
f = cec13_func(x, 26) - 1200;
if f < 1e-8
	f = 0;
end
end

