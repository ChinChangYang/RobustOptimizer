mainfile	= 'cec14D50_jadebin_Q4_201404170333';
subfile		= 'cec14D50_jade_s_Q4_201404170541';
mainfilename = sprintf('%s.mat', mainfile);
xlsfilename = sprintf('%s.xlsx', mainfile);
subfilename = sprintf('%s.mat', subfile);

% Generate Measurements
load(mainfilename);
close all;
allfvals(allfvals <= 1e-8) = 0;
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
for i = 1 : nfuncs
	allfvalssorted(:, :, i) = allfvals(:, sortindices(:, :, i), i);
end
errmedian	= allfvalssorted(:, round(0.5 * (end + 1)), :);
errmedian	= reshape(errmedian, nprogress, nfuncs)';
[NP, ~]		= size(allout{1, 1}.FC);
q			= zeros(nruns, NP, nprogress, nfuncs);
for i = 1 : nfuncs
	for j = 1 : nruns
		q(j, :, :, i) = allout{j, i}.FC;
	end
end
qsorted = q;
for i = 1 : nfuncs
	qsorted(:, :, :, i) = q(sortindices(:, :, i), :, :, i);
end
qmedian		= qsorted(round(0.5 * (end + 1)), :, :, :);
qmedianmean = mean(qmedian, 2);
qmedianmean = reshape(qmedianmean, nprogress, nfuncs)';
qmedianmax	= max(qmedian, [], 2);
qmedianmax	= reshape(qmedianmax, nprogress, nfuncs)';
qmedian		= reshape(qmedian, NP, nprogress, nfuncs);
qmediansum	= sum(qmedian > solverOptions.Q);
qmediansum = reshape(qmediansum, nprogress, nfuncs)';
fes			= allout{1, 1}.fes;
G			= allout{1, 1}.G;

if isfield(allout{1, 1}, 'MF')
	MF		= zeros(nruns, NP, nprogress, nfuncs);
	for i = 1 : nfuncs
		for j = 1 : nruns
			MF(j, :, :, i) = allout{j, i}.MF;
		end
	end
	MF		= reshape(MF, nruns * NP, nprogress, nfuncs);	
end

if isfield(allout{1, 1}, 'MCR')
	MCR		= zeros(nruns, NP, nprogress, nfuncs);
	for i = 1 : nfuncs
		for j = 1 : nruns
			MCR(j, :, :, i) = allout{j, i}.MCR;
		end
	end
	MCR		= reshape(MCR, nruns * NP, nprogress, nfuncs);	
end

% Save to Excel
xlswrite(xlsfilename, {'Mean'}, 'Error', 'A1');
xlswrite(xlsfilename, errmean, 'Error', 'A2:A31');
xlswrite(xlsfilename, {'St. D.'}, 'Error', 'B1');
xlswrite(xlsfilename, errstd, 'Error', 'B2:B31');
xlswrite(xlsfilename, {'SR'}, 'Error', 'C1');
xlswrite(xlsfilename, succrate, 'Error', 'C2:C31');
xlswrite(xlsfilename, G, 'Convergence', 'A1:U1');
xlswrite(xlsfilename, errmedian, 'Convergence', 'A2:U31');
xlswrite(xlsfilename, G, 'q > Q', 'A1:U1');
xlswrite(xlsfilename, qmediansum, 'q > Q', 'A2:U31');

% Convergence Graph (Example)
figure;
semilogy(G, errmedian(18, :), 'k');
title('Function 18');
xlabel('Generation');
ylabel('Solution Error');

% Dynamic of q value (Example)
figure;
boxplot(qmedian(:, :, 18), G, 'colors', 'k', 'plotstyle','compact');
title('Function 18');
xlabel('Generation');
ylabel('q');

% Dynamic of the number of q > Q
figure;
plot(G, qmediansum(18, :), 'k');
title('Function 18');
xlabel('Generation');
ylabel('Number of q > Q');

% Dynamic of MF value (Example)
if isfield(allout{1, 1}, 'MF')
	figure;
	boxplot(MF(:, :, 18), G, 'colors', 'k', 'plotstyle','compact');
	title('Function 18');
	xlabel('Generation');
	ylabel('MF');
end

% Dynamic of MCR value (Example)
if isfield(allout{1, 1}, 'MCR')
	figure;
	boxplot(MCR(:, :, 18), G, 'colors', 'k', 'plotstyle','compact');
	title('Function 18');
	xlabel('Generation');
	ylabel('MCR');
end

% Wilcoxon Rank Sum Test
load(mainfilename);
allfvals(allfvals <= 1e-8) = 0;
A = reshape(allfvals(end, :, :), nruns, nfuncs);
load(subfilename);
allfvals(allfvals <= 1e-8) = 0;
B = reshape(allfvals(end, :, :), nruns, nfuncs);
w			= ranksumtest(A, B);
POSITIVE	= sum(w=='+');
EQUAL		= sum(w=='=');
NEGATIVE	= sum(w=='-');