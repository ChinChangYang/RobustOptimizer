function [c, ceq] = sainz_c1(x, y)
%SAINZ_C1 Constraint of Problem 1 (Sainz et al., 2008)

nx = numel(x);
twice_nx = nx + nx;
c = zeros(twice_nx, 1);
c(1 : nx) = y - x .* (x + 6.28);
c(nx + 1 : twice_nx) = y - x .* (x - 6.28);
ceq = 0;
end
