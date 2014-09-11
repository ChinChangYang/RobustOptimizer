function [c, ceq] = fminmax_c1(x, y)
c = x(1) - x(2);
ceq = abs(y(1) - y(2));
end
