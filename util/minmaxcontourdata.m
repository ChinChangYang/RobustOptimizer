function [XX, YY, ZZ] = minmaxcontourdata(D, lb1, ub1, lb2, ub2, fitfun, options)
if nargin <= 6
	options = [];
end

defaultOptions.Samples_One_Dimension = 100;
defaultOptions.Samples_Multi_Dimension = 10;
defaultOptions.ScaleType = 'sqrt';
options = setdefoptions(options, defaultOptions);
n1 = options.Samples_One_Dimension;
n2 = options.Samples_Multi_Dimension;
scaleType = options.ScaleType;
	
if D == 1
	my = 0.5 * (lb2 + ub2);
	vx = linspace(lb1(1), ub1(1), n1);
	vy = linspace(lb2(1), ub2(1), n1);
	[XX, YY] = meshgrid(vx, vy);
	ZZ = zeros(numel(vx), numel(vy));
	for i = 1 : numel(vx)
		parfor j = 1 : numel(vy)
			y = my;
			y(1) = YY(i, j);
			ZZ(i, j) = feval(fitfun, XX(i, j), y);
		end
	end
else
	mx = 0.5 * (lb1 + ub1);
	vx = linspace(lb1(1), ub1(1), n2);
	vy = linspace(lb1(2), ub1(2), n2);
	[XX, YY] = meshgrid(vx, vy);
	ZZ = zeros(numel(vx), numel(vy));
	for i = 1 : numel(vx)
		parfor j = 1 : numel(vy)
			x = mx;
			x(1:2) = [XX(i, j); YY(i, j)];
			fitfunX2i = @(y) -feval(fitfun, x, y);
			XX2 = jadebin(fitfunX2i, lb2, ub2, numel(lb2)^2 * 1e2);
			ZZ(i, j) = feval(fitfun, x, XX2);
		end
	end
end

if strcmp(scaleType, 'sqrt')
	ZZ = sqrt(ZZ - min(ZZ(:)) + 1);
elseif strcmp(scaleType, 'log')	
	ZZ = -log(-ZZ - min(-ZZ(:)) + 1);
elseif strcmp(scaleType, 'sqrtsqrt')
	ZZ = (ZZ - min(ZZ(:)) + 1).^0.25;
end
end
