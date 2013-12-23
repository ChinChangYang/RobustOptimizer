function f = prtb_rast(x,y)
%PRTB_RAST Perturbed Rastrigin's function
f = -(10 + sum((y - 0.5 * x).^2) - sum(10 * x * cos(2 * pi *(y - 0.5*x))));

end
