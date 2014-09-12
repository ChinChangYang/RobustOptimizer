function testminmaxsolver
%TESTMINMAXSOLVER Test a minmax solver
if matlabpool('size') == 0
	matlabpool('open');
end

startTime = tic;
close all;
rng(1, 'twister');
% solver = 'minmaxdegl';
% solver = 'mmjade_pce';
% solver = 'mmdeb1b_pce';
% solver = 'mmder1b_pce';
solver = 'mmshade_pce';
% solver = 'mmsade_pce';
% fitfun = 'sainz_f1';
% nonlcon = 'sainz_c1';
% fitfun = 'sainz_f2';
% nonlcon = 'sainz_c2';
% fitfun = 'lu_f1';
% nonlcon = 'lu_c1';
fitfun = 'fminmax_f2';
nonlcon = 'fminmax_c2';
D_Min = 2;
D_Max = 2;
D = D_Min + D_Max;
maxfunevals = D * 5e6;
solverOptions1.NP = 30;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.9;
solverOptions1.Display = 'iter';
solverOptions1.RecordPoint = 1000;
% solverOptions1.InnerSolver = 'shade';
solverOptions1.nonlcon = nonlcon;
solverOptions1.innerMaxIter = 200;
solverOptions1.migrateFactor = 0.7;
solverOptions1.ConstraintHandling = 'EpsilonMethod';
solverOptions1.EpsilonValue = 1e-6;
solverOptions1.EarlyStop = 'auto';
% solverOptions1.InnerSolver = innerSolver;
solverOptions2.NP = 30;
solverOptions2.F = 0.7;
solverOptions2.CR = 0.9;
solverOptions2.TolStagnationIteration = 20;
solverOptions2.ConstraintHandling = 'EpsilonMethod';
solverOptions2.EpsilonValue = 1e-6;
solverOptions2.EarlyStop = 'auto';
% lb1 = -3.14 * ones(D_Min, 1);
% ub1 = 3.14 * ones(D_Min, 1);
% lb2 = -3.14 * ones(D_Max, 1);
% ub2 = 3.14 * ones(D_Max, 1);
% lb1 = 0 * ones(D_Min, 1);
% ub1 = 6 * ones(D_Min, 1);
% lb2 = 2 * ones(D_Max, 1);
% ub2 = 8 * ones(D_Max, 1);
% lb1 = -5 * ones(D_Min, 1);
% ub1 = 5 * ones(D_Min, 1);
% lb2 = -5 * ones(D_Max, 1);
% ub2 = 5 * ones(D_Max, 1);
% lb1 = 0 * ones(D_Min, 1);
% ub1 = 10 * ones(D_Min, 1);
% lb2 = 0 * ones(D_Max, 1);
% ub2 = 10 * ones(D_Max, 1);
% lb1 = -100 * ones(D_Min, 1);
% ub1 = 100 * ones(D_Min, 1);
% lb2 = -100 * ones(D_Max, 1);
% ub2 = 100 * ones(D_Max, 1);
lb1 = -2 * ones(D_Min, 1);
ub1 = 2 * ones(D_Min, 1);
lb2 = -4 * ones(D_Max, 1);
ub2 = 4 * ones(D_Max, 1);
% lb1 = -2e50 * ones(D_Min, 1);
% ub1 = 2e50 * ones(D_Min, 1);
% lb2 = -4e50 * ones(D_Max, 1);
% ub2 = 4e50 * ones(D_Max, 1);
[xminmax1, xminmax2, fminmax, out] = ...
    feval(solver, fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, solverOptions1, solverOptions2);
% fprintf('FEs = \n');
% disp(out.fes(end));
fprintf('xminmax1 = \n');
disp(xminmax1);
fprintf('xminmax2 = \n');
disp(xminmax2);
fprintf('fminmax = %.4E\n', fminmax);
fprintf('xminmax1 - 0.1 = %.4E\n', norm(xminmax1 - 0.1));
fprintf('countcon = %.4E\n', out.countcon);

optima = feval(fitfun);
[~, n_optima] = size(optima);
distance = zeros(n_optima, 1);

for i_optima = 1 : n_optima
	distance(i_optima) = norm([xminmax1; xminmax2] - optima(:, i_optima));
end
			
[min_distance, ~] = min(distance);			
fprintf('xminmax1 - xoptim = %.4E\n', min_distance);
% fprintf('out.bestever.fminmax = \n');
% disp(out.bestever.fminmax);
% fprintf('out.bestever.xminmax1 = \n');
% disp(out.bestever.xminmax1);
% fprintf('out.bestever.xminmax2 = \n');
% disp(out.bestever.xminmax2);
% fprintf('xmean = \n');
% disp(out.xmean(:, end));
figure;
hold off;
plot(out.fes, out.fmean);
xlabel('FEs');
title('function values');
figure;
semilogy(out.fes, out.fstd);
xlabel('FEs');
title('Std. of function values');
figure;
hold off;
for i = 1 : numel(out.xmean(:, 1))
    plot(out.fes, out.xmean(i, :), getlinespec(i));
    hold on;
end
xlabel('FEs');
title('Mean of X solutions');
figure;
hold off;
for i = 1 : numel(out.xstd(:, 1))
	semilogy(out.fes, out.xstd(i, :), getlinespec(i));
    hold on;
end
xlabel('FEs');
title('St. D. of X solutions');
if isfield(out, 'distancemin')
	figure;
	hold off;
	semilogy(out.fes, out.distancemin, 'r');
	hold on;
	semilogy(out.fes, out.distancemax, 'r');
	semilogy(out.fes, out.distancemean, 'b');
	semilogy(out.fes, out.distancemedian, 'g');
	xlabel('FEs');
	title('Distances between each pair of X');
end

if isfield(out, 'cond');
	figure;
	hold off;
	semilogy(out.fes, out.cond);
	xlabel('FEs');
	title('Condition number');
end

if isfield(out, 'innerFstd');
	figure;
	hold off;
	semilogy(out.fes, out.innerFstd);
	xlabel('FEs');
	title('Std. f in inner states');
end

if isfield(out, 'innerMeanXstd')
	figure;
	hold off;
	semilogy(out.fes, out.innerMeanXstd);
	xlabel('FEs');
	title('Mean of std. X in inner states');
end

if isfield(out, 'angle')
	figure;
	NP = solverOptions1.NP;
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
end

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

if isfield(out, 'alpha')
	figure;
	plot(out.fes, out.alpha);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('alpha');
end

if isfield(out, 'successRate')
	figure;
	plot(out.fes, out.successRate);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('Success Rate');
end

if isfield(out, 'X_Converged_FEs')
	figure;
	plot(out.fes, out.X_Converged_FEs);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('X_Converged_FEs');
end

if isfield(out, 'U_Converged_FEs')
	figure;
	plot(out.fes, out.U_Converged_FEs);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('U_Converged_FEs');
end
toc(startTime);
filename = sprintf('data_%s_%s.mat', solver, fitfun);
save(filename, 'xminmax1', 'xminmax2', 'fminmax', 'out');
end
