function [c, ceq] = liu_c1(x, y)
%LIU_C1 Constraint of Problem 4 (Liu et al., 1998)

c = zeros(3, 1);
c(1) = x(1).^2 + x(2).^2 - 100;
c(2) = y(1) - x(1);
c(3) = y(2) - x(2);
ceq = 0;
end
