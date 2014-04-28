clear;
load('filenames_201404241706.mat');
nfilenames = numel(filenames);
nmainalgo = nfilenames / 2;
filenames_all = cell(2, nmainalgo);
for iSolver = 1 : 2 : nfilenames
	filenames_all{1, 0.5 * (iSolver + 1)} = filenames{iSolver};
end
for iSolver = 2 : 2 : nfilenames
	filenames_all{2, 0.5 * iSolver} = filenames{iSolver};
end
D = 30;
Q = 46;
tablefilename		= sprintf('CEC14_D%d_Q%d_TABLE.xlsx', D, Q);
convfilename		= sprintf('CEC14_D%d_Q%d_CONV.xlsx', D, Q);
qdynfliename		= sprintf('CEC14_D%d_Q%d_QDYN.xlsx', D, Q);
xstdfilename		= sprintf('CEC14_D%d_Q%d_XSTD.xlsx', D, Q);
solver_all			= cell(2, nmainalgo);
ranksumtestsheet	= 'Rank Sum Test';

xlswrite(tablefilename, ...
	{'Method 1', 'Method 2', 'POSITIVE', 'EQUAL', 'NEGATIVE'}, ...
	ranksumtestsheet, ...
	'A1:E1');

load(filenames_all{1, 1});
xlswrite(tablefilename, ...
	{'NP'; 'F'; 'CR'; 'Q'; 'D'; 'RUNS'; 'MaxFEs'}, ...
	'Param', ...
	'A1:A7');
xlswrite(tablefilename, ...
	{solverOptions.NP; ...
	solverOptions.F; ...
	solverOptions.CR; ...
	solverOptions.Q; ...
	measureOptions.Dimension; ...
	measureOptions.Runs; ...
	measureOptions.MaxFunEvals}, ...
	'Param', ...
	'B1:B7');

for iSolver = 1 : nmainalgo
	for iSV = 1 : 2
		load(filenames_all{iSV, iSolver});
		solver_all{iSV, iSolver} = solver;
		
		% Generate Measurements
		allfvals(allfvals <= 1e-8) = 0; %#ok<*SAGROW>
		errmean		= mean(allfvals(end, :, :), 2);
		errmean		= errmean(:);
		errstd		= std(allfvals(end, :, :), [], 2);
		errstd		= errstd(:);
		succrate	= mean(allfvals(end, :, :) <= 1e-8);
		succrate	= succrate(:);
		compcomplex	= (T2 - T1) / T0;
		[nprogress, nruns, nfuncs] = size(allfvals);
		[~, sortindices] = sort(allfvals(end, :, :), 2);
		allfvalssorted = allfvals;
		for j = 1 : nfuncs
			allfvalssorted(:, :, j) = allfvals(:, sortindices(:, :, j), j);
		end
		errmedian	= allfvalssorted(:, round(0.5 * (end + 1)), :);
		errmedian	= reshape(errmedian, nprogress, nfuncs)';
		[NP, ~]		= size(allout{1, 1}.FC);
		q			= zeros(nruns, NP, nprogress, nfuncs);
		for j = 1 : nfuncs
			for k = 1 : nruns
				q(k, :, :, j) = allout{k, j}.FC;
			end
		end
		qsorted = q;
		for j = 1 : nfuncs
			qsorted(:, :, :, j) = q(sortindices(:, :, j), :, :, j);
		end
		qmedian		= qsorted(round(0.5 * (end + 1)), :, :, :);
		qmedianmean = mean(qmedian, 2);
		qmedianmean = reshape(qmedianmean, nprogress, nfuncs)';
		qmedianmax	= max(qmedian, [], 2);
		qmedianmax	= reshape(qmedianmax, nprogress, nfuncs)';
		qmedian		= reshape(qmedian, NP, nprogress, nfuncs);
		qmediansum	= sum(qmedian > solverOptions.Q);
		qmediansum	= reshape(qmediansum, nprogress, nfuncs)';
		xstd		= zeros(nruns, nprogress, nfuncs);
		for j = 1 : nfuncs
			for k = 1 : nruns
				xstd(k, :, j) = mean(allout{k, j}.xstd);
			end
		end
		xstdsorted = xstd;
		for j = 1 : nfuncs
			xstdsorted(:, :, j) = xstd(sortindices(:, :, j), :, j);
		end
		xstdmedian	= xstdsorted(round(0.5 * (end + 1)), :, :);
		xstdmedian	= reshape(xstdmedian, nprogress, nfuncs)';
		
		% Save to Excel
		fes			= allout{1, 1}.fes;
		G			= allout{1, 1}.G;
		switch iSV
			case 1
				xlswrite(tablefilename, {solver}, solver, 'A1');
				xlswrite(tablefilename, {'Mean'}, solver, 'A2');
				xlswrite(tablefilename, errmean, solver, 'A3:A32');
				xlswrite(tablefilename, {'St. D.'}, solver, 'B2');
				xlswrite(tablefilename, errstd, solver, 'B3:B32');
				xlswrite(tablefilename, {'SR'}, solver, 'C2');
				xlswrite(tablefilename, succrate, solver, 'C3:C32');
			case 2
				solver_wo_sv = solver_all{1, iSolver};
				xlswrite(tablefilename, {solver}, solver_wo_sv, 'D1');
				xlswrite(tablefilename, {'Mean'}, solver_wo_sv, 'D2');
				xlswrite(tablefilename, errmean, solver_wo_sv, 'D3:D32');
				xlswrite(tablefilename, {'St. D.'}, solver_wo_sv, 'E2');
				xlswrite(tablefilename, errstd, solver_wo_sv, 'E3:E32');
				xlswrite(tablefilename, {'SR'}, solver_wo_sv, 'F2');
				xlswrite(tablefilename, succrate, solver_wo_sv, 'F3:F32');
		end
		
		xlswrite(convfilename, G, solver, 'A1:U1');
		xlswrite(convfilename, errmedian, solver, 'A2:U31');
		xlswrite(qdynfliename, G, solver, 'A1:U1');
		xlswrite(qdynfliename, qmediansum, solver, 'A2:U31');
		xlswrite(xstdfilename, G, solver, 'A1:U1');
		xlswrite(xstdfilename, xstdmedian, solver, 'A2:U31');
	end
		
	% Wilcoxon Rank Sum Test
	load(filenames_all{1, iSolver});
	allfvals(allfvals <= 1e-8) = 0;
	A = reshape(allfvals(end, :, :), nruns, nfuncs);
	load(filenames_all{2, iSolver});
	allfvals(allfvals <= 1e-8) = 0;
	B = reshape(allfvals(end, :, :), nruns, nfuncs);
	w			= ranksumtest(A, B);
	POSITIVE	= sum(w=='+');
	EQUAL		= sum(w=='=');
	NEGATIVE	= sum(w=='-');
	
	solver_wo_sv = solver_all{1, iSolver};
	solver_w_sv = solver_all{2, iSolver};
	data = {solver_wo_sv, solver_w_sv, POSITIVE, EQUAL, NEGATIVE};
	
	range = sprintf('A%d:E%d', iSolver + 1, iSolver + 1);
	xlswrite(tablefilename, ...
		data, ...
		ranksumtestsheet, ...
		range);
	
	for i = 1 : numel(w)
		range = sprintf('G%d', i + 2);
		xlswrite(tablefilename, w(i), solver_wo_sv, range);
	end
end
