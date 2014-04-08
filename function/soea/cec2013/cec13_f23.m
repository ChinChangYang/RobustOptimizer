function f = cec13_f23(x)
% CEC13_F23 CEC'2013 function 23
% Composition Function 3
% n = 3
% g1: Rotated Schwefel's Function f15
% g2: Rotated Schwefel's Function f15
% g3: Rotated Schwefel's Function f15
x = reshape(x, numel(x), 1);
f = cec13_func(x, 23) - 900;
end

