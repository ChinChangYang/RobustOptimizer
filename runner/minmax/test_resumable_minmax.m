function [norm_resum, norm_no_resum] = test_resumable_minmax
if matlabpool('size') == 0
    matlabpool('open');
end

startTime = tic;
nRuns = 40;
T = 50;
D1 = 1;
D2 = 1;
DF1 = 10;
DF2 = 10;

norm_resum = zeros(1, nRuns);
norm_no_resum = zeros(1, nRuns);
solver = 'minmaxtcjadebin';
fitfun = 'fminmax_f2';
maxfunevals = D1 * D2 * DF1 * DF2 * 5e2;
solverOptions1.dimensionFactor = DF1;
solverOptions1.RecordPoint = 0;
solverOptions2.dimensionFactor = DF2;
solverOptions2.RecordPoint = 0;
lb1 = -2e60 * ones(D1, 1);
ub1 = 2e60 * ones(D1, 1);
lb2 = -4e60 * ones(D2, 1);
ub2 = 4e60 * ones(D2, 1);

parfor iRuns = 1 : nRuns
	%% Resume
	rng(iRuns, 'twister');
	fes = 0;
	
	[xminmax1, ~, ~, out] = ...
		feval(solver, fitfun, 1/T * maxfunevals, lb1, ub1, lb2, ub2, ...
		solverOptions1, solverOptions2);
	
	fes = fes + out.fes(end);
	
	while fes < maxfunevals && out.fes(end) > 0
		resumeOptions1 = solverOptions1;
		resumeOptions1.initial = out.final;
		[xminmax1, ~, ~, out] = ...
			feval(solver, fitfun, 1/T * maxfunevals, lb1, ub1, lb2, ub2, ...
		resumeOptions1, solverOptions2);
	
		fes = fes + out.fes(end);
	end
	
	norm_resum(iRuns) = norm(xminmax1 - feval(fitfun));
	
	%% No Resume
	rng(iRuns, 'twister');
	[xminmax1, ~, ~, ~] = ...
		feval(solver, fitfun, fes, lb1, ub1, lb2, ub2, ...
		solverOptions1, solverOptions2);
	
	norm_no_resum(iRuns) = norm(xminmax1 - feval(fitfun));
	
	%% Report
	fprintf('Run: %d, Done.\n', iRuns);
end

norm_resum = sort(norm_resum);
fprintf('Mean of norm_resum: %.4E\n', mean(norm_resum));
fprintf('St. D of norm_resum: %.4E\n', std(norm_resum));
fprintf('Min of norm_resum: %.4E\n', norm_resum(1));
fprintf('Q1 of norm_resum: %.4E\n', norm_resum(round(0.25 * nRuns)));
fprintf('Q2 of norm_resum: %.4E\n', norm_resum(round(0.5 * nRuns)));
fprintf('Q3 of norm_resum: %.4E\n', norm_resum(round(0.75 * nRuns)));
fprintf('Max of norm_resum: %.4E\n', norm_resum(end));

norm_no_resum = sort(norm_no_resum);
fprintf('Mean of norm_no_resum: %.4E\n', mean(norm_no_resum));
fprintf('St. D of norm_no_resum: %.4E\n', std(norm_no_resum));
fprintf('Min of norm_no_resum: %.4E\n', norm_no_resum(1));
fprintf('Q1 of norm_no_resum: %.4E\n', norm_no_resum(round(0.25 * nRuns)));
fprintf('Q2 of norm_no_resum: %.4E\n', norm_no_resum(round(0.5 * nRuns)));
fprintf('Q3 of norm_no_resum: %.4E\n', norm_no_resum(round(0.75 * nRuns)));
fprintf('Max of norm_no_resum: %.4E\n', norm_no_resum(end));

[~, h, ~] = ranksum(norm_resum, norm_no_resum);

if h == 1
	fprintf('norm_resum is significantly difference with norm_no_resum\n');
else
	fprintf('No significant difference between norm_resum and norm_no_resum\n');
end

toc(startTime);
end

