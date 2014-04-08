function f = cec13_f22(x)
% CEC13_F22 CEC'2013 function 22
% Composition Function 2
% n = 3
% g1: Schwefel's Function f14
% g2: Schwefel's Function f14
% g3: Schwefel's Function f14
x = reshape(x, numel(x), 1);
f = cec13_func(x, 22) - 800;
end

