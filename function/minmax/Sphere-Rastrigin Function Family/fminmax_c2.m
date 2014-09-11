function [c, ceq] = fminmax_c2(x, y)
c = 0;
ceq = [abs(x(1) - x(2)); abs(y(1) - y(2))];
end
