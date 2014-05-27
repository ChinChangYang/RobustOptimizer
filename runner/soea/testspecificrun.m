%TESTSPECIFICRUN Test a specific experiment of certain solver, test
%function, maximal function evaluations.
startTime = tic;
clear;
close all;
rng('default');
% load('InitialX.mat');
load('InitialX_CEC11.mat');
solver = 'jade_sps';
fitfun = 'cec11_f12';
D = 26;
maxfunevals = D * 1e4;
solverOptions.nonlcon = [];
solverOptions.Q = 64;
solverOptions.NP = 5 * D;
solverOptions.ftarget = 1e-8;
solverOptions.Restart = 0;
solverOptions.Display = 'iter';
solverOptions.RecordPoint = 21;
solverOptions.Noise = false;
% solverOptions.initial.X = eval(sprintf('XD%dNP%d', ...
% 	D, ...
% 	solverOptions.NP));
solverOptions.initial.X = eval(sprintf('XF%dNP%d', ...
	3, ...
	solverOptions.NP));
% lb = -5 * ones(D, 1);
% ub = 5 * ones(D, 1);
% lb = -100 * ones(D, 1);
% ub = 100 * ones(D, 1);
% lb = -6.4 * ones(D, 1);
% ub = 6.35 * ones(D, 1);
% lb = zeros(20, 1);
% ub = 4 * pi * ones(20, 1);
[lb, ub] = getlimit_messenger;
[xmin, fmin, out] = ...
    feval(solver, fitfun, lb, ub, maxfunevals, solverOptions);
% fprintf('out.bestever.xmin = \n');
% disp(out.bestever.xmin);
% if isfield(out, 'xmean')
% 	fprintf('xmean = \n');
% 	disp(out.xmean(:, end));
% end
% fprintf('out.bestever.fmin = \n');
% disp(out.bestever.fmin);
% fprintf('xmin = \n');
% disp(xmin);
% fprintf('fes = %.4E\n', out.fes(end));
fprintf('fmin = %.4E\n', fmin);
% if strncmp('bbob12', fitfun, 6)
% 	fprintf('fmin - fopt = %.4E\n', fmin - feval(fitfun, 'xopt'));
% end
if isfield(out, 'stopflag')
	fprintf('stopflag = %s\n', out.stopflag);
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
	loglog(out.G, out.distancemean);
else
	semilogy(out.G, out.distancemean);
end
xlabel('Generation');
title('Mean of distances between the target vectors and the centroid');
subplot(236);
hold off;
semilogy(out.fes, out.cond);
xlabel('FEs');
title('Condition number');

if isfield(solverOptions, 'dimensionFactor')
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
end

if isfield(out, 'mu_F')
	figure;
	plot(out.fes, out.mu_CR);
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

if isfield(out, 'NP')
	figure;
	semilogy(out.fes, out.NP);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('NP');
end

if isfield(out, 'converg_rate')
	figure;
	semilogy(out.fes, out.converg_rate);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('converg rate');
end

if isfield(out, 'm')
	figure;
	plot(out.fes, out.m);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('m');
end

if isfield(out, 'geomean_Sconvrate')
	figure;
	plot(out.fes, out.geomean_Sconvrate);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('geomean_Sconvrate');
end

if isfield(out, 'mu_H')
	figure;
	semilogy(out.fes, out.mu_H);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_H');
end

if isfield(out, 'MF')
	figure;
	if isvector(out.MF)
		plot(out.fes, out.MF);
		title(sprintf('Solve %s by %s', fitfun, solver));
		xlabel('FEs');
		ylabel('MF');
	else
		boxplot(out.MF, out.G, 'colors', 'k', 'plotstyle','compact');
		title(sprintf('Solve %s by %s', fitfun, solver));
		xlabel('Generation');
		ylabel('MF');
	end
end

if isfield(out, 'MCR')
	figure;
	if isvector(out.MCR)
		plot(out.fes, out.MCR);
		title(sprintf('Solve %s by %s', fitfun, solver));
		xlabel('FEs');
		ylabel('MCR');
	else
		boxplot(out.MCR, out.G, 'colors', 'k', 'plotstyle','compact');
		title(sprintf('Solve %s by %s', fitfun, solver));
		xlabel('Generation');
		ylabel('MCR');
	end
end

if isfield(out, 'MR')
	figure;
	plot(out.fes, out.MR);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('Mean of R');
end

if isfield(out, 'countStagnation')
	figure;
	plot(out.fes, out.countStagnation);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('countStagnation');
end

if isfield(out, 'MG')
	figure;
	plot(out.fes, out.MG);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('MG');
end

if isfield(out, 'MPMAX')
	figure;
	plot(out.fes, out.MPMAX);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('MPMAX');
end

if isfield(out, 'FCMEDIAN')
	figure;
	hold on;
	plot(out.fes, out.FC1Q, 'g');
	plot(out.fes, out.FCMEDIAN, 'r');
	plot(out.fes, out.FC3Q, 'b');
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('Recent Consecutive Unsuccessful Trial Vectors');
	legend('1Q', 'MEDIAN', '3Q');
% 	print(sprintf('%s.tiff', fitfun), '-dtiff');
end

if isfield(out, 'FCMEAN')
	figure;
	errorbar(out.fes, out.FCMEAN, out.FCSTD);
	axis([0, out.fes(end), 0, max(out.FCMEAN) + 2 * max(out.FCSTD)]);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('Function Evaluations');
	ylabel('Recent Consecutive Unsuccessful Trial Vectors');
% 	print(sprintf('%s.tiff', fitfun), '-dtiff');
end

if isfield(out, 'FC')
	figure;
	boxplot(out.FC, out.G, 'colors', 'k', 'plotstyle','compact');
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('Generation');
	ylabel('Recent Consecutive Unsuccessful Trial Vectors');
% 	print(sprintf('%s.tiff', fitfun), '-dtiff');
end

if isfield(out, 'FC')
	figure;
	plot(out.G, mean(out.FC));
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('Generation');
	ylabel('Recent Consecutive Unsuccessful Trial Vectors');
% 	print(sprintf('%s.tiff', fitfun), '-dtiff');
end

figure;
semilogy(out.fes, mean(out.xstd));
title(sprintf('Solve %s by %s', fitfun, solver));
xlabel('FEs');
ylabel('St. D. of Target Vectors');

if isfield(out, 'S_FC') && isfield(out, 'S_mFC')
	fprintf('mean(out.S_FC) = %f\n', mean(out.S_FC));
	fprintf('mean(out.S_mFC) = %f\n', mean(out.S_mFC));
	if mean(out.S_FC) < mean(out.S_mFC)
		fprintf('mean(out.S_FC) < mean(out.S_mFC)\n');
	else
		fprintf('mean(out.S_FC) >= mean(out.S_mFC)\n');
	end
end

if isfield(out, 'mFC') && isfield(out, 'mSFC')
	plot(out.G, out.mFC, 'b');
	hold on;
	plot(out.G, out.mSFC, 'r');
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('Generation');
	ylabel('Number');
end

% toc(startTime);