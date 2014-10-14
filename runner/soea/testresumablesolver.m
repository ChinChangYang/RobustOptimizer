function testresumablesolver
if matlabpool('size') == 0
    matlabpool('open');
end

startTime = tic;
nRuns = 40;
T = 50;

fmin_resum = zeros(1, nRuns);
fmin_no_resum = zeros(1, nRuns);
solver = 'cmaes';
fitfun = 'bbob12_f1';
D = 10;
maxfunevals = D * 4e3;
solverOptions.dimensionFactor = 5;
solverOptions.NP = 5 * D;
solverOptions.TolX = 100 * eps;
solverOptions.TolFun = eps;
solverOptions.Display = 'off';
solverOptions.RecordPoint = 0;
lb = -5e50 * ones(D, 1);
ub = 5e50 * ones(D, 1);

parfor iRuns = 1 : nRuns
	% Solver with resuming
	rng(iRuns, 'twister');
	resumeOptions = solverOptions;
	resumeOptions.usefunevals = 1 / T * maxfunevals;
	
	[~, ~, out] = ...
		feval(solver, fitfun, lb, ub, maxfunevals, resumeOptions);
	
	for t = 1 : T - 1	
		resumeOptions.initial = out.final;
		resumeOptions.usefunevals = (t + 1) / T * maxfunevals;
		[~, fmin_resum(iRuns), out] = ...
			feval(solver, fitfun, lb, ub, maxfunevals, resumeOptions);
	end
		
	% Solver without resuming
	rng(iRuns, 'twister');
	[~, fmin_no_resum(iRuns), ~] = ...
		feval(solver, fitfun, lb, ub, out.fes(end), solverOptions);
	
	fprintf('Run: %d, Done.\n', iRuns);
end

fmin_resum = sort(fmin_resum);
fprintf('Mean of fmin_resum: %.4E\n', mean(fmin_resum));
fprintf('St. D of fmin_resum: %.4E\n', std(fmin_resum));
fprintf('Min of fmin_resum: %.4E\n', fmin_resum(1));
fprintf('Q1 of fmin_resum: %.4E\n', fmin_resum(round(0.25 * nRuns)));
fprintf('Q2 of fmin_resum: %.4E\n', fmin_resum(round(0.5 * nRuns)));
fprintf('Q3 of fmin_resum: %.4E\n', fmin_resum(round(0.75 * nRuns)));
fprintf('Max of fmin_resum: %.4E\n', fmin_resum(end));

fmin_no_resum = sort(fmin_no_resum);
fprintf('Mean of fmin_no_resum: %.4E\n', mean(fmin_no_resum));
fprintf('St. D of fmin_no_resum: %.4E\n', std(fmin_no_resum));
fprintf('Min of fmin_no_resum: %.4E\n', fmin_no_resum(1));
fprintf('Q1 of fmin_no_resum: %.4E\n', fmin_no_resum(round(0.25 * nRuns)));
fprintf('Q2 of fmin_no_resum: %.4E\n', fmin_no_resum(round(0.5 * nRuns)));
fprintf('Q3 of fmin_no_resum: %.4E\n', fmin_no_resum(round(0.75 * nRuns)));
fprintf('Max of fmin_no_resum: %.4E\n', fmin_no_resum(end));

[~, h, ~] = ranksum(fmin_resum, fmin_no_resum);

if h == 1
	fprintf('fmin_resume is significantly difference with fmin_no_resume\n');
else
	fprintf('No significant difference between fmin_resume and fmin_no_resume\n');
end

toc(startTime);
end
