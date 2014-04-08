function f = cec13_f14(x)
% CEC13_F14 CEC'2013 function 14
% Schwefel's Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 14) + 100;
end

