function [XX, YY, ZZ] = preparecontourdata(D, lb, ub, fitfun)
% Initialize contour data
if D >= 2
	vx = linspace(lb(1), ub(1), 50);
	vy = linspace(lb(2), ub(2), 50);
	[XX, YY] = meshgrid(vx, vy);
	ZZ = zeros(numel(vx), numel(vy));
	v = 0.5 * (lb + ub);
	for i = 1 : numel(vx)
		for j = 1 : numel(vy)
			v(1:2, 1) = [XX(i, j); YY(i, j)];
			ZZ(i, j) = feval(fitfun, v);
		end
	end
	ZZ = log(ZZ - min(ZZ(:)) + 1);
else
	XX = linspace(lb, ub, 2000);
	YY = zeros(1, numel(XX));
	for i = 1 : numel(XX)
		YY(i) = feval(fitfun, XX(i));
	end
	ZZ = [];
end
end

