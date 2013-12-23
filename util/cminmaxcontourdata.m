function [XX, YY, ZZ, CC] = cminmaxcontourdata(D, lb1, ub1, lb2, ub2, fitfun, options)
%CONSMINMAXCONTOURDATA Constraint min-max contour data

if nargin <= 6
	options = [];
end

[XX, YY, ZZ] = minmaxcontourdata(D, lb1, ub1, lb2, ub2, fitfun, options);

defaultOptions.nonlcon = [];
defaultOptions.Samples_One_Dimension = 100;
options = setdefoptions(options, defaultOptions);
nonlcon = options.nonlcon;
n1 = options.Samples_One_Dimension;

if D > 1 || isempty(nonlcon)
	CC = [];
else
	my = 0.5 * (lb2 + ub2);
	vx = linspace(lb1(1), ub1(1), n1);
	vy = linspace(lb2(1), ub2(1), n1);
	[XX, YY] = meshgrid(vx, vy);
	CC = zeros(numel(vx), numel(vy));	
	for i = 1 : numel(vx)
		parfor j = 1 : numel(vy)
			y = my;
			y(1) = YY(i, j);
			c = feval(nonlcon, XX(i, j), y);
			CC(i, j) = sum(c(c > 0));
		end
	end
	
	CC(CC < 0) = -1;
	CC(CC > 0) = 1;
end
end

