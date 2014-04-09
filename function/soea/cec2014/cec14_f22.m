function f = cec14_f22(x)
% CEC14_F22 CEC'2014 function 22
% Hybrid Function 6
% N = 5
% p = [0.1, 0.2, 0.2, 0.2, 0.3]
% g1 : Katsuura Function f10
% g2 : HappyCat Function f11
% g3: Expanded Griewank・s plus Rosenbrock・s Function f13
% g4: Modified Schwefel・s Function f9
% g5: Ackley・s Function f5
x = reshape(x, numel(x), 1);
f = cec14_func(x, 22) - 2200;
end
