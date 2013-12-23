function f = prtb_quad(x, y)
% PTRB_QUAD Perturbed quadratic function
f = sum(x.^2) - 2 * sum((y-x).^2);
end
