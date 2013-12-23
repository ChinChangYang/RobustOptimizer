function f = cec13_f27(x)
% CEC13_F27 CEC'2013 function 27
% Composition Function 7
% n = 5
% g1: Rotated Griewank's Function f10
% g2: Rotated Rastrigins Function f12
% g3: Rotated Schwefel's Function f15
% g4: Rotated Weierstrass Function f9
% g5: Sphere Function f1
x = reshape(x, numel(x), 1);
f = cec13_func(x, 27) - 1300;
if f < 1e-8
	f = 0;
end
end

