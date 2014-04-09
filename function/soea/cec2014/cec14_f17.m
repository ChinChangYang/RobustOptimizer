function f = cec14_f17(x)
% CEC14_F17 CEC'2014 function 17
% Hybrid Function 1
% p = [0.3,0.3,0.4]
% g1 : Modified Schwefel's Function f9
% g2 : Rastrigin¡¦s Function f8
% g3: High Conditioned Elliptic Function f1
x = reshape(x, numel(x), 1);
f = cec14_func(x, 17) - 1700;
end
