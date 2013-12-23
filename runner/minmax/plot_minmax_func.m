function plot_minmax_func(fitfun)
tic;

if nargin == 0
	fitfun = 'fminmax_f8';
end

rng('default');
close all;
lb = [-2; -4];
ub = [2; 4];

% X = (lb+ub)/2
sample = 1000;
vy = linspace(lb(2), ub(2), sample);
f = zeros(1, numel(vy));
optimum = feval(fitfun);
x_optim = optimum;

for i = 1 : numel(vy)
	f(i) = feval(fitfun, x_optim, vy(i));
end

figure(1);
plot(vy, f, 'k');
title(sprintf('x = %.1f', x_optim));
xlabel('y');
ylabel('f');
pause(0.01);
print(sprintf('%s_1', fitfun), '-dtiff');

% Y = arg(y) maxminmax f(x, y)
sample = 50;
vx = linspace(lb(1), ub(1), sample);
f = zeros(1, numel(vx));

for i = 1 : numel(vx)
	fitfunXi = @(y) -feval(fitfun, vx(i), y);
	dfmaxfunevals = 5e3;
	defaultDF = 20;
	options.dimensionFactor = defaultDF;
	[~, f(i), out] = jadebin(fitfunXi, lb(2), ub(2), dfmaxfunevals, options);
	
	factor = 1;
	while max(out.xstd(:, end)) > 1e-2 && factor < 100
		factor = 2 * factor;
		maxfunevals = dfmaxfunevals * factor;
		options.dimensionFactor = round(defaultDF * factor ^ 0.25);
		[~, f(i), out] = jadebin(fitfunXi, lb(2), ub(2), maxfunevals, options);
	end
	
	if factor > 1
		fprintf('Max converged at factor: %d\n', factor);
	end
	
	f(i) = -f(i);
end

figure(2);
plot(vx, f, 'k');
title('y = arg(y) maxminmax f(x, y)');
xlabel('x');
ylabel('f');
pause(0.01);
print(sprintf('%s_2', fitfun), '-dtiff');

% X-Y contour
scaleType1 = 'none';
contourline = 50;
sample = 200;
vx = linspace(lb(1), ub(1), sample);
vy = linspace(lb(2), ub(2), sample);
[XX, YY] = meshgrid(vx, vy);
FF = zeros(numel(vx), numel(vy));
for i = 1 : numel(vx)
	for j = 1 : numel(vy)
		FF(i, j) = feval(fitfun, XX(i, j), YY(i, j));
	end
end

if strcmp(scaleType1, 'sqrt')
	FF = sqrt(FF - min(FF(:)) + 1);
elseif strcmp(scaleType1, 'log')
	FF = -log(-FF - min(-FF(:)) + 1);
end

figure(3);
contour(XX, YY, FF, contourline);
cmap = colormap('gray');
cmap = cmap(end:-1:1, :);
colormap(cmap);
colorbar;
xlabel('x');
ylabel('y');
pause(0.01);
print(sprintf('%s_3', fitfun), '-dtiff');
toc;
end
