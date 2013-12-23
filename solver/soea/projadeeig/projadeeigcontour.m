function projadeeigcontour
%PROJADEEIGCONTOUR Contour plot of Pro JADE/eig
load projadeeigcontour;
contour(X, Y, Z, 30, 'LineWidth', 2);
xlabel('dimensionFactor');
ylabel('delta_F');
title('Contour plot of solution errors of Pro JADE/eig for D = 2, M = 100');
end
