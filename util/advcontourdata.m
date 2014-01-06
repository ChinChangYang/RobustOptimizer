function [XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun)
% Initialize contour data
if D >= 3
	N = 50;
	cD = min(3, D - 1);
	XX = zeros(N, N, cD);
	YY = zeros(N, N, cD);
	ZZ = zeros(N, N, cD);
	
	for d = 1 : cD
		v = 0.5 * (lb + ub);
		vx = linspace(lb(d), ub(d), N);
		vy = linspace(lb(d + 1), ub(d + 1), N);
		[XX(:, :, d), YY(:, :, d)] = meshgrid(vx, vy);
		for i = 1 : N
			for j = 1 : N
				v([d, d+1], 1) = [XX(i, j, d); YY(i, j, d)];
				ZZ(i, j, d) = feval(fitfun, v);
			end
		end
		ZZ(:, :, d) = log(ZZ(:, :, d) - min(min(ZZ(:, :, d))) + 1);
	end
elseif D == 2
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

