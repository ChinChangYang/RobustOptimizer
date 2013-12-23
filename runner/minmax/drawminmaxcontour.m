function drawminmaxcontour
addprojectpath;
close all;
fitfun = 'minmax_func6';
lb = [-100; -100];
ub = [100; 100];
vx = linspace(lb(1), ub(1), 50);
vy = linspace(lb(2), ub(2), 50);
[XX, YY] = meshgrid(vx, vy);
ZZ = zeros(numel(vx), numel(vy));
n = 50;
scaleType1 = 'sqrt';
scaleType2 = 'sqrt';

%% Draw 1-D minimizer and 1-D maximizer contour
for i = 1 : numel(vx)
	for j = 1 : numel(vy)
		ZZ(i, j) = feval(fitfun, XX(i, j), YY(i, j));
	end
end

if strcmp(scaleType1, 'sqrt')
	ZZ = sqrt(ZZ - min(ZZ(:)) + 1);
elseif strcmp(scaleType1, 'log')
	ZZ = -log(-ZZ - min(-ZZ(:)) + 1);
end

figure(1);
contour(XX, YY, ZZ, n);
colorbar;
xlabel('x');
ylabel('y');
pause(0.01);
print(sprintf('%s_1', fitfun), '-dtiff');

%% Draw 2-D minimizer and 2-D maximizer contour
for i = 1 : numel(vx)
	for j = 1 : numel(vy)
		ZZ(i, j) = feval(fitfun, lb, [XX(i, j); YY(i, j)]);
	end
end
if strcmp(scaleType1, 'sqrt')
	ZZ = sqrt(ZZ - min(ZZ(:)) + 1);
elseif strcmp(scaleType1, 'log')
	ZZ = -log(-ZZ - min(-ZZ(:)) + 1);
end
figure(2);
contour(XX, YY, ZZ, n);
colorbar;
title('X=lb');
xlabel('y_1');
ylabel('y_2');
pause(0.01);
print(sprintf('%s_2', fitfun), '-dtiff');

for i = 1 : numel(vx)
	for j = 1 : numel(vy)
		ZZ(i, j) = feval(fitfun, 0.75 * lb + 0.25 * ub, [XX(i, j); YY(i, j)]);
	end
end
if strcmp(scaleType1, 'sqrt')
	ZZ = sqrt(ZZ - min(ZZ(:)) + 1);
elseif strcmp(scaleType1, 'log')
	ZZ = -log(-ZZ - min(-ZZ(:)) + 1);
end
figure(3);
contour(XX, YY, ZZ, n);
colorbar;
title('X=0.75*lb+0.25*ub');
xlabel('y_1');
ylabel('y_2');
pause(0.01);
print(sprintf('%s_3', fitfun), '-dtiff');

for i = 1 : numel(vx)
	for j = 1 : numel(vy)
		ZZ(i, j) = feval(fitfun, 0.5 * (lb + ub), [XX(i, j); YY(i, j)]);
	end
end
if strcmp(scaleType1, 'sqrt')
	ZZ = sqrt(ZZ - min(ZZ(:)) + 1);
elseif strcmp(scaleType1, 'log')
	ZZ = -log(-ZZ - min(-ZZ(:)) + 1);
end
figure(4);
contour(XX, YY, ZZ, n);
colorbar;
title('X=0.5*lb+0.5*ub');
xlabel('y_1');
ylabel('y_2');
pause(0.01);
print(sprintf('%s_4', fitfun), '-dtiff');

for i = 1 : numel(vx)
	for j = 1 : numel(vy)
		ZZ(i, j) = feval(fitfun, 0.25 * lb + 0.75 * ub, [XX(i, j); YY(i, j)]);
	end
end
if strcmp(scaleType1, 'sqrt')
	ZZ = sqrt(ZZ - min(ZZ(:)) + 1);
elseif strcmp(scaleType1, 'log')
	ZZ = -log(-ZZ - min(-ZZ(:)) + 1);
end
figure(5);
contour(XX, YY, ZZ, n);
colorbar;
title('X=0.25*lb+0.75*ub');
xlabel('y_1');
ylabel('y_2');
pause(0.01);
print(sprintf('%s_5', fitfun), '-dtiff');

for i = 1 : numel(vx)
	for j = 1 : numel(vy)
		ZZ(i, j) = feval(fitfun, ub, [XX(i, j); YY(i, j)]);
	end
end
if strcmp(scaleType1, 'sqrt')
	ZZ = sqrt(ZZ - min(ZZ(:)) + 1);
elseif strcmp(scaleType1, 'log')
	ZZ = -log(-ZZ - min(-ZZ(:)) + 1);
end
figure(6);
contour(XX, YY, ZZ, n);
colorbar;
title('X=ub');
xlabel('y_1');
ylabel('y_2');
pause(0.01);
print(sprintf('%s_6', fitfun), '-dtiff');

options.Samples_Multi_Dimension = 30;
options.ScaleType = scaleType2;
[XX, YY, ZZ] = minmaxcontourdata(2, lb, ub, lb, ub, fitfun, options);
figure(7);
contour(XX, YY, ZZ, n);
title('Y=arg max_y f(x,y)');
xlabel('x_1');
ylabel('x_2');
pause(0.01);
print(sprintf('%s_7', fitfun), '-dtiff');
end
