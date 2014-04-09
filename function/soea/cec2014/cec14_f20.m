function f = cec14_f20(x)
% CEC14_F20 CEC'2014 function 20
% Hybrid Function 4
% N = 4
% p = [0.2, 0.2, 0.3, 0.3]
% g1 : HGBat Function f12
% g2 : Discus Function f3
% g3: Expanded Griewank¡¦s plus Rosenbrock¡¦s Function f13
% g4: Rastrigin¡¦s Function f8
x = reshape(x, numel(x), 1);
f = cec14_func(x, 20) - 2000;
end
