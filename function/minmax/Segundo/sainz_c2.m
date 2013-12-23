function [c, ceq] = sainz_c2(x, y)
%SAINZ_C2 Constraint of Problem 2 (Sainz et al., 2008)

nx = numel(x);
twice_nx = nx + nx;
c = zeros(twice_nx, 1);
c(1 : nx)		= sum(-(x - 5).^2 - (y - 3).^2 + 4);
c(nx + 1 : twice_nx) = sum((x - 5).^2 + (y - 3).^2 - 16);
ceq = 0;
end
