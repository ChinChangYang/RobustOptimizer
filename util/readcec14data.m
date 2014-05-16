function readcec14data(metafilename, metafilename_sps)
load(metafilename);
filenames_o = filenames;
load(metafilename_sps);
filenames_sps = filenames;
[nQ, nA] = size(filenames_sps);
pn = zeros(nA, nQ);
Q = zeros(1, nQ);
for iQ = 1 : nQ
	load(filenames_sps{iQ, 1});
	D					= measureOptions.Dimension;
	Q(iQ)				= solverOptions.Q;
	tablefilename		= sprintf('CEC14_D%d_Q%d_TABLE.xlsx', D, Q(iQ));
	convfilename		= sprintf('CEC14_D%d_Q%d_CONV.xlsx', D, Q(iQ));
	qdynfilename		= sprintf('CEC14_D%d_Q%d_QDYN.xlsx', D, Q(iQ));
	contfilename		= sprintf('CEC14_D%d_Q%d_CONT.xlsx', D, Q(iQ));
	qbarfilename		= sprintf('CEC14_D%d_Q%d_QBAR.xlsx', D, Q(iQ));
	pnfilename			= sprintf('CEC14_D%d_PN.xlsx', D);
	solver_all			= cell(2, nA);
	ranksumtestsheet	= 'Rank Sum Test';
	
	xlswrite(tablefilename, ...
		{'Method 1', 'Method 2', 'POSITIVE', 'EQUAL', 'NEGATIVE', 'PN'}, ...
		ranksumtestsheet, ...
		'A1:F1');
	
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
	
	for iA = 1 : nA
		[solver, errmean, errstd, succrate, compcomplex, errmedian, ...
			qmediansum, distancemedian, fes, G, qmedianmean] = ...
			readonedata(filenames_o{iA});
		
		solver_all{1, iA} = solver;
		
		xlswrite(tablefilename, {solver}, solver, 'A1');
		xlswrite(tablefilename, {'Mean'}, solver, 'A2');
		xlswrite(tablefilename, errmean, solver, 'A3:A32');
		xlswrite(tablefilename, {'St. D.'}, solver, 'B2');
		xlswrite(tablefilename, errstd, solver, 'B3:B32');
		xlswrite(tablefilename, {'SR'}, solver, 'C2');
		xlswrite(tablefilename, succrate, solver, 'C3:C32');
		
		
		xlswrite(convfilename, G, solver, 'A1:U1');
		xlswrite(convfilename, errmedian, solver, 'A2:U31');
		xlswrite(qdynfilename, G, solver, 'A1:U1');
		xlswrite(qdynfilename, qmediansum, solver, 'A2:U31');
		xlswrite(contfilename, G, solver, 'A1:U1');
		xlswrite(contfilename, distancemedian, solver, 'A2:U31');
		xlswrite(qbarfilename, G, solver, 'A1:U1');
		xlswrite(qbarfilename, qmedianmean, solver, 'A2:U31');
		
		[solver, errmean, errstd, succrate, compcomplex, errmedian, ...
			qmediansum, distancemedian, fes, G, qmedianmean] = ...
			readonedata(filenames_sps{iQ, iA});
		
		solver_all{2, iA} = solver;
		
		solver_wo_sv = solver_all{1, iA};
		xlswrite(tablefilename, {solver}, solver_wo_sv, 'D1');
		xlswrite(tablefilename, {'Mean'}, solver_wo_sv, 'D2');
		xlswrite(tablefilename, errmean, solver_wo_sv, 'D3:D32');
		xlswrite(tablefilename, {'St. D.'}, solver_wo_sv, 'E2');
		xlswrite(tablefilename, errstd, solver_wo_sv, 'E3:E32');
		xlswrite(tablefilename, {'SR'}, solver_wo_sv, 'F2');
		xlswrite(tablefilename, succrate, solver_wo_sv, 'F3:F32');
		
		xlswrite(convfilename, G, solver, 'A1:U1');
		xlswrite(convfilename, errmedian, solver, 'A2:U31');
		xlswrite(qdynfilename, G, solver, 'A1:U1');
		xlswrite(qdynfilename, qmediansum, solver, 'A2:U31');
		xlswrite(contfilename, G, solver, 'A1:U1');
		xlswrite(contfilename, distancemedian, solver, 'A2:U31');
		xlswrite(qbarfilename, G, solver, 'A1:U1');
		xlswrite(qbarfilename, qmedianmean, solver, 'A2:U31');
		
		% Wilcoxon Rank Sum Test		
		load(filenames_o{iA});
		[~, nruns_A, nfuncs_A] = size(allfvals);
		load(filenames_sps{iQ, iA});
		[~, nruns_B, nfuncs_B] = size(allfvals);
		nruns = min(nruns_A, nruns_B);
		
		% Check nfuncs
		if nfuncs_A ~= nfuncs_B
			error('nfuncs_A ~= nfuncs_B');
		else
			nfuncs = nfuncs_A;
		end
		
		load(filenames_o{iA});
		allfvals(allfvals <= 1e-8) = 0; %#ok<*SAGROW>
		A = reshape(allfvals(end, :, :), nruns_A, nfuncs);
		load(filenames_sps{iQ, iA});
		[~, nruns_B, nfuncs] = size(allfvals);
		allfvals(allfvals <= 1e-8) = 0;
		B = reshape(allfvals(end, :, :), nruns_B, nfuncs);
		w			= ranksumtest(A(1:nruns, :), B(1:nruns, :));
		POSITIVE	= sum(w=='+');
		EQUAL		= sum(w=='=');
		NEGATIVE	= sum(w=='-');
		pn(iA, iQ)	= POSITIVE - NEGATIVE;
		
		solver_wo_sv = solver_all{1, iA};
		solver_w_sv = solver_all{2, iA};
		data = {solver_wo_sv, solver_w_sv, POSITIVE, EQUAL, NEGATIVE, pn(iA, iQ)};
		
		range = sprintf('A%d:F%d', iA + 1, iA + 1);
		xlswrite(tablefilename, ...
			data, ...
			ranksumtestsheet, ...
			range);
		
		range = 'G3:G32';
		xlswrite(tablefilename, w, solver_wo_sv, range);
		
		% Overall P-N
		xlswrite(pnfilename, ...
			[Q(iQ), POSITIVE, EQUAL, NEGATIVE, pn(iA, iQ)], ...
			solver, ...
			sprintf('A%d', iQ + 1));
	end
end

for iA = 1 : nA	
	load(filenames_sps{1, iA});
	xlswrite(pnfilename, ...
		{'Q', 'POSITIVE', 'EQUAL', 'NEGATIVE', 'P-N'}, ...
		solver, ...
		'A1');
	
	xlswrite(pnfilename, ...
		{solver}, ...
		'ALL', ...
		sprintf('A%d', iA + 1));
	
	xlswrite(pnfilename, ...
		pn(iA, :), ...
		'ALL', ...
		sprintf('B%d', iA + 1));
end

xlswrite(pnfilename, ...
	{'Solver'}, ...
	'ALL', ...
	'A1');

xlswrite(pnfilename, ...
	Q, ...
	'ALL', ...
	'B1');

xlswrite(pnfilename, ...
	{'SUM'}, ...
	'ALL', ...
	sprintf('A%d', nA + 2));

xlswrite(pnfilename, ...
	sum(pn), ...
	'ALL', ...
	sprintf('B%d', nA + 2));
end

function [solver, errmean, errstd, succrate, compcomplex, errmedian, ...
	qmediansum, distancemedian, fes, G, qmedianmean] = readonedata(filename)

load(filename);

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
distance	= zeros(nruns, nprogress, nfuncs);
for j = 1 : nfuncs
	for k = 1 : nruns
		distance(k, :, j) = allout{k, j}.distancemean;
	end
end
distancesorted = distance;
for j = 1 : nfuncs
	distancesorted(:, :, j) = distance(sortindices(:, :, j), :, j);
end
distancemedian	= distancesorted(round(0.5 * (end + 1)), :, :);
distancemedian	= reshape(distancemedian, nprogress, nfuncs)';

fes			= allout{1, 1}.fes;
G			= allout{1, 1}.G;
end