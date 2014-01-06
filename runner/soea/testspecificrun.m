function testspecificrun
%TESTSPECIFICRUN Test a specific experiment of certain solver, test
%function, maximal function evaluations.
startTime = tic;
close all;
solver = 'jadebin';
fitfun = 'cec13_f2';
D = 5;
maxfunevals = D * 1e4;
solverOptions.nonlcon = [];
solverOptions.dimensionFactor = 5;
solverOptions.F = 0.7;
solverOptions.CR = 0.5;
solverOptions.R = 0.5;
solverOptions.TolX = 1e3 * eps;
solverOptions.TolFun = 1e3 * eps;
solverOptions.ftarget = -Inf;
solverOptions.Restart = 0;
solverOptions.Display = 'iter';
solverOptions.RecordPoint = 1000;
solverOptions.Noise = false;
% lb = -5e50 * ones(D, 1);
% ub = 5e50 * ones(D, 1);
lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);
% lb = -6.4 * ones(D, 1);
% ub = 6.35 * ones(D, 1);
% lb = zeros(20, 1);
% ub = 4 * pi * ones(20, 1);
% [lb, ub] = getlimit_messenger;
[xmin, fmin, out] = ...
    feval(solver, fitfun, lb, ub, maxfunevals, solverOptions);
fprintf('out.bestever.xmin = \n');
disp(out.bestever.xmin);
if isfield(out, 'xmean')
	fprintf('xmean = \n');
	disp(out.xmean(:, end));
end
fprintf('out.bestever.fmin = \n');
disp(out.bestever.fmin);
fprintf('xmin = \n');
disp(xmin);
fprintf('fmin = %.4E\n', fmin);
if strncmp('bbob12', fitfun, 6)
	fprintf('fmin - fopt = %.4E\n', fmin - feval(fitfun, 'xopt'));
end
figure(2);
subplot(231);
hold off;
lowerBoundF = min([out.fmean, out.fmin]);
semilogy(out.fes, out.fmean - lowerBoundF, 'b');
hold on;
semilogy(out.fes, out.fmin - lowerBoundF, 'r');
xlabel('FEs');
title('function values');
subplot(232);
if solverOptions.Noise
	loglog(out.fes, out.fstd);
else
	semilogy(out.fes, out.fstd);
end
xlabel('FEs');
title('Std. of function values');
subplot(233);
hold off;
for i = 1 : numel(out.xmean(:, 1))
    plot(out.fes, out.xmean(i, :), getlinespec(i));
    hold on;
end
xlabel('FEs');
title('Mean of X solutions');
subplot(234);
hold off;
for i = 1 : numel(out.xstd(:, 1))
	if solverOptions.Noise
		loglog(out.fes, out.xstd(i, :), getlinespec(i));
	else
		semilogy(out.fes, out.xstd(i, :), getlinespec(i));
	end
    hold on;
end
xlabel('FEs');
title('Std. of X solutions');
subplot(235);
hold off;
if solverOptions.Noise
	loglog(out.fes, out.distancemin, 'r');
	hold on;
	loglog(out.fes, out.distancemax, 'r');
	loglog(out.fes, out.distancemean, 'b');
	loglog(out.fes, out.distancemedian, 'g');
else
	semilogy(out.fes, out.distancemin, 'r');
	hold on;
	semilogy(out.fes, out.distancemax, 'r');
	semilogy(out.fes, out.distancemean, 'b');
	semilogy(out.fes, out.distancemedian, 'g');
end
xlabel('FEs');
title('Distances between each pair of X');
subplot(236);
hold off;
semilogy(out.fes, out.cond);
xlabel('FEs');
title('Condition number');

figure;
NP = D * solverOptions.dimensionFactor;
Gmax = maxfunevals / NP;
alternative_angle = out.angle;
alternative_angle(out.angle > pi/4) = out.angle(out.angle > pi/4) - pi/2;
if std(out.angle) < std(alternative_angle)
	semilogx(out.fes / NP, out.angle / pi * 180, 'k');
	axis([0, Gmax, 0, 90]);
else
	semilogx(out.fes / NP, alternative_angle / pi * 180, 'k');
	axis([0, Gmax, -45, 45]);
end
xlabel('Generation');
ylabel('Angle (degree)');
title('Angle between the natural basis and the eigenvector basis');

if isfield(out, 'mu_F')
	figure;
	semilogy(out.fes, out.mu_F);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_F');
end

if isfield(out, 'mu_CR')
	figure;
	plot(out.fes, out.mu_CR);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_CR');
end

if isfield(out, 'mu_R')
	figure;
	plot(out.fes, out.mu_R);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_R');
end

if isfield(out, 'mu_G')
	figure;
	plot(out.fes, out.mu_G);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_G');
end
toc(startTime);
end
