startTime = tic;
clear;
close all;
load('InitialX.mat');
solver = 'debest1bin_e';
fitfuns = {'cec14_f2', 'cec14_f4'};
D = 30;
maxfunevals = D * 1e4;
solverOptions.Q = 21;
solverOptions.NP = 5 * D;
solverOptions.ftarget = 1e-8;
solverOptions.Restart = 0;
solverOptions.Display = 'off';
solverOptions.RecordPoint = 21;
solverOptions.Noise = false;
solverOptions.initial.X = eval(sprintf('XD%dNP%d', ...
	D, ...
	solverOptions.NP));

lb = -100 * ones(D, 1);
ub = 100 * ones(D, 1);
linespecs = {'k', 'k-.'; ...
	'k--', 'k:'};
h1 = figure;
h2 = figure;

for iFitfun = 1 : numel(fitfuns)	
	fitfun = fitfuns{iFitfun};
	rng('default');
	
	[xmin, fmin, out] = ...
		feval(solver, fitfun, lb, ub, maxfunevals, solverOptions);
	
	figure(h1);
	semilogy(out.G, out.distancemean, linespecs{iFitfun});
	hold on;
		
	figure(h2);
	plot(out.G, out.mFC, linespecs{iFitfun, 1});
	hold on;
	plot(out.G, out.mSFC, linespecs{iFitfun, 2});
end

figure(h1);
xlabel('\fontname{Times New Roman}Generation');
ylabel('\fontname{Times New Roman}Distance');
title({'\fontname{Times New Roman}Mean of distances', ...
	'\fontname{Times New Roman}between target vectors and centroid'});
legend('\fontname{Times New Roman}\itf_{\rm2}', ...
	'\fontname{Times New Roman}\itf_{\rm4}');

figure(h2);
% legend({'\fontname{Times New Roman}\itmq_G \rm(\itf_{\rm2}\rm)', ...
% 	'\fontname{Times New Roman}\itmsq_G \rm(\itf_{\rm2}\rm)', ...
% 	'\fontname{Times New Roman}\itmq_G \rm(\itf_{\rm4}\rm)', ...
% 	'\fontname{Times New Roman}\itmsq_G \rm(\itf_{\rm4}\rm)'});

title({'\fontname{Times New Roman}Number of recent', ...
	'\fontname{Times New Roman}consecutive unsuccessful trial vectors'});
xlabel('\fontname{Times New Roman}Generation');
ylabel('\fontname{Times New Roman}Number');
legend({'$$\bar{q}_G\,\,\,\,\,\,\,\,\,(f_2)$$', ...
	'$$\bar{q}_G^{(succ)}\,(f_2)$$', ...
	'$$\bar{q}_G\,\,\,\,\,\,\,\,\,(f_4)$$', ...
	'$$\bar{q}_G^{(succ)}\,(f_4)$$'}, ...
	'Interpreter', 'latex');