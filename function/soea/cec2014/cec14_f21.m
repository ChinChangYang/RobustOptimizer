function f = cec14_f21(x)
% CEC14_F21 CEC'2014 function 21
% Hybrid Function 5
% N = 5
% p = [0.1, 0.2, 0.2, 0.2, 0.3]
% g1 : Scaffer¡¦s F6 Function:f14
% g2 : HGBat Function f12
% g3: Rosenbrock¡¦s Function f4
% g4: Modified Schwefel¡¦s Function f9
% g5: High Conditioned Elliptic Function f1
x = reshape(x, numel(x), 1);
f = cec14_func(x, 21) - 2100;
end
