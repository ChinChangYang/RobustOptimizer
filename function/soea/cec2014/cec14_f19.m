function f = cec14_f19(x)
% CEC14_F19 CEC'2014 function 19
% Hybrid Function 3
% N = 4
% p = [ 0.2, 0.2, 0.3, 0.3]
% g1 : Griewank¡¦s Function f7
% g2 : Weierstrass Function f6
% g3: Rosenbrock¡¦s Function f4
% g4: Scaffer¡¦s F6 Function:f14
x = reshape(x, numel(x), 1);
f = cec14_func(x, 19) - 1900;
end
