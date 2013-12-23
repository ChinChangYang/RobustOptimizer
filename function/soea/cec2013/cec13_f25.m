function f = cec13_f25(x)
% CEC13_F25 CEC'2013 function 25
% Composition Function 5
% n = 3
% g1: Rotated Schwefel's Function f15
% g2: Rotated Rastrigin's Function f12
% g3: Rotated Weierstrass Function f9
%
% All settings are same as Composition Function 4, except
% sigma = [10, 30, 50]
x = reshape(x, numel(x), 1);
f = cec13_func(x, 25) - 1100;
if f < 1e-8
	f = 0;
end
end

