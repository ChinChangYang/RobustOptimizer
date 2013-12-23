function ZETA_Param_Analysis
if matlabpool('size') == 0
	matlabpool('open');
end

startTime = tic;
close all;
rng(1, 'twister');
solver = 'mmmdeb1b_c';
fitfun = 'maxminmax_f1';
D1 = 1;
D2 = 1;
D3 = 1;
% D = D1 + D2 + D3;
maxfunevals = 1e8;
solverOptions1.dimensionFactor = 10;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.9;
solverOptions1.TolX = eps;
solverOptions1.TolFun = 0;
solverOptions1.Display = 'off';
solverOptions1.RecordPoint = 1;
solverOptions2.dimensionFactor = 10;
solverOptions2.F = 0.7;
solverOptions2.CR = 0.9;
solverOptions3.dimensionFactor = 10;
solverOptions3.F = 0.7;
solverOptions3.CR = 0.9;
lb1 = -1e50 * ones(D1, 1);
ub1 = 1e50 * ones(D1, 1);
lb2 = -2e50 * ones(D2, 1);
ub2 = 2e50 * ones(D2, 1);
lb3 = -4e50 * ones(D3, 1);
ub3 = 4e50 * ones(D3, 1);

nRuns = 8;
zeta = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

for i = 1 : numel(zeta)
	err = zeros(nRuns, 1);
	fes = zeros(nRuns, 1);
	solverOptions1.zeta = zeta(i);
	
	parfor j = 1 : nRuns
		[xbest1, ~, ~, ~, out] = ...
			feval(solver, fitfun, ...
			maxfunevals, lb1, ub1, lb2, ub2, lb3, ub3, solverOptions1, ...
			solverOptions2, solverOptions3);
		
		err(j) = norm(xbest1 - 0.1);
		fes(j) = out.fes(end);
	end
	
	err_mean = mean(err);
	err_std = std(err);
	fes_mean = mean(fes);
	fes_std = std(fes);
	
	fprintf('err_mean = \n');
	disp(err_mean);
	fprintf('err_std = \n');
	disp(err_std);	
	fprintf('fes_mean = \n');
	disp(fes_mean);
	fprintf('fes_std = \n');
	disp(fes_std);
	
	save(sprintf('err_maxminmax_%s_%s.mat', solver, datestr(now, 'yyyymmddHHMM')), ...
		'err', 'err_mean', 'err_std', ...
		'fes', 'fes_mean', 'fes_std', ...
		'solverOptions1', 'solverOptions2', 'solverOptions3');
	fprintf('zeta: %s\n', zeta(i));
	fprintf('done\n');
end

toc(startTime);
end

