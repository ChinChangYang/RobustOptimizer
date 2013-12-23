function [c, ceq] = lu_c1(x, y)
%LU_C1 Constraint of Problem 3 (Lu et al., 2008)

nx = numel(x);
c = zeros(nx, 1);
c(1 : nx) = sum(-x.^2 - y.^2 + 25);
ceq = 0;
end
