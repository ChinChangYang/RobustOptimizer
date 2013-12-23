function plot_maxminmax_func(fitfun)
tic;

if nargin == 0
	fitfun = 'maxminmax_f8';
end

rng('default');
close all;
lb = [-1; -2; -4];
ub = [1; 2; 4];
contourline = 30;
scaleType1 = 'none';
scaleType2 = 'none';

% X = (lb+ub)/2; Y = (lb+ub)/2
sample = 1000;
vz = linspace(lb(3), ub(3), sample);
f = zeros(1, numel(vz));
optimum = feval(fitfun);
x_optim = optimum;
y_optim = optimum;

for i = 1 : numel(vz)
	f(i) = feval(fitfun, x_optim, y_optim, vz(i));
end

figure(1);
plot(vz, f, 'k');
title(sprintf('(x, y) = (%.1f, %.1f)', x_optim, y_optim));
xlabel('z');
ylabel('f');
pause(0.01);
print(sprintf('%s_1', fitfun), '-dtiff');

% X = 0.5 * lb + 0.5 * ub
sample = 100;
vy = linspace(lb(2), ub(2), sample);
vz = linspace(lb(3), ub(3), sample);
[YY, ZZ] = meshgrid(vy, vz);
FF = zeros(numel(vy), numel(vz));

for i = 1 : numel(vy)
	for j = 1 : numel(vz)
		FF(i, j) = feval(fitfun, x_optim, YY(i, j), ZZ(i, j));
	end
end

if strcmp(scaleType1, 'sqrt')
	FF = sqrt(FF - min(FF(:)) + 1);
elseif strcmp(scaleType1, 'log')
	FF = -log(-FF - min(-FF(:)) + 1);
end

figure(2);
contour(YY, ZZ, FF, contourline);
cmap = colormap('gray');
cmap = cmap(end:-1:1, :);
colormap(cmap);
colorbar;
title(sprintf('x = %.1f', x_optim));
xlabel('y');
ylabel('z');
pause(0.01);
print(sprintf('%s_2', fitfun), '-dtiff');

% Z = arg(z) maxminmax f(x,y,z)
sample = 30;
contourline = 30;
vx = linspace(lb(1), ub(1), sample);
vy = linspace(lb(2), ub(2), sample);
[XX, YY] = meshgrid(vx, vy);
FF = zeros(numel(vx), numel(vy));

for i = 1 : numel(vx)
	for j = 1 : numel(vy)		
		fitfunXiYj = @(z) -feval(fitfun, XX(i, j), YY(i, j), z);
		dfmaxfunevals = 5e3;
		defaultDF = 20;
		options.dimensionFactor = defaultDF;
		[~, FF(i, j), out] = jadebin(fitfunXiYj, lb(3), ub(3), dfmaxfunevals, options);	
		
		factor = 1;
		while max(out.xstd(:, end)) > 1e-2 && factor < 100
			factor = 2 * factor;
			maxfunevals = dfmaxfunevals * factor;
			options.dimensionFactor = round(defaultDF * factor ^ 0.25);
			[~, FF(i, j), out] = jadebin(fitfunXiYj, lb(3), ub(3), maxfunevals, options);
		end
		
		if factor > 1
			fprintf('Max converged at factor: %d\n', factor);
		end
		
		FF(i, j) = -FF(i, j);	
	end
end

if strcmp(scaleType2, 'sqrt')
	FF = sqrt(FF - min(FF(:)) + 1);
elseif strcmp(scaleType2, 'log')
	FF = -log(-FF - min(-FF(:)) + 1);
end

figure(3);
contour(XX, YY, FF, contourline);
cmap = colormap('gray');
cmap = cmap(end:-1:1, :);
colormap(cmap);
colorbar;
title('z = arg(z) maxminmax f(x, y, z)');
xlabel('x');
ylabel('y');
pause(0.01);
print(sprintf('%s_3', fitfun), '-dtiff');

% [Y, Z] = arg(y,z) maxminmax f(x,y,z)
sample = 30;
vx = linspace(lb(1), ub(1), sample);
f = zeros(1, numel(vx));

for i = 1 : numel(vx)	
	fitfunXi = @(y, z) feval(fitfun, vx(i), y, z);
	dfmaxfunevals = 2e5;
	defaultDF = 15;
	options1.dimensionFactor = defaultDF;
	options2.dimensionFactor = defaultDF;
	[~, ~, f(i), out] = minmaxtcjadebin(fitfunXi, dfmaxfunevals, ...
		lb(2), ub(2), lb(3), ub(3), options1, options2);
	
	factor = 1;
	while max(out.xstd(:, end)) > 1e-2 && factor < 100
		factor = 2 * factor;
		options1.dimensionFactor = round(defaultDF * factor ^ 0.25);
		options2.dimensionFactor = round(defaultDF * factor ^ 0.25);
		[~, ~, f(i), out] = minmaxtcjadebin(fitfunXi, round(dfmaxfunevals * factor), ...
			lb(2), ub(2), lb(3), ub(3), options1, options2);
	end
	
	if factor > 1
		fprintf('Minmax converged at factor: %d\n', factor);
	end
end

figure(4);
plot(vx, f, 'k');
title('(y, z) = arg(y,z) maxminmax f(x, y, z)');
xlabel('x');
ylabel('f');
pause(0.01);
print(sprintf('%s_4', fitfun), '-dtiff');

% Y = arg(y) maxminmax f(x,y,z)
sample = 30;
contourline = 30;
vx = linspace(lb(1), ub(1), sample);
vz = linspace(lb(3), ub(3), sample);
[XX, ZZ] = meshgrid(vx, vz);
FF = zeros(numel(vx), numel(vz));

for i = 1 : numel(vx)
	for j = 1 : numel(vz)		
		fitfunXiZj = @(y) feval(fitfun, XX(i, j), y, ZZ(i, j));
		dfmaxfunevals = 5e3;
		defaultDF = 20;
		options.dimensionFactor = defaultDF;
		[~, FF(i, j), out] = jadebin(fitfunXiZj, lb(2), ub(2), dfmaxfunevals, options);	
		
		factor = 1;
		while max(out.xstd(:, end)) > 1e-2 && factor < 100
			factor = 2 * factor;
			maxfunevals = dfmaxfunevals * factor;
			options.dimensionFactor = round(defaultDF * factor ^ 0.25);
			[~, FF(i, j), out] = jadebin(fitfunXiZj, lb(2), ub(2), maxfunevals, options);
		end
		
		if factor > 1
			fprintf('Min converged at factor: %d\n', factor);
		end	
	end
end

if strcmp(scaleType1, 'sqrt')
	FF = sqrt(FF - min(FF(:)) + 1);
elseif strcmp(scaleType1, 'log')
	FF = -log(-FF - min(-FF(:)) + 1);
end

figure(5);
contour(XX, ZZ, FF, contourline);
cmap = colormap('gray');
cmap = cmap(end:-1:1, :);
colormap(cmap);
colorbar;
title('y = arg(y) maxminmax f(x, y, z)');
xlabel('x');
ylabel('z');
pause(0.01);
print(sprintf('%s_5', fitfun), '-dtiff');
toc;
end
