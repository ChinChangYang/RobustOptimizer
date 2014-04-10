mainfile	= 'complete_cec13_201404101001';
subfile		= 'complete_cec13_201404101001';
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
q			= reshape(q, nruns * NP, nprogress, nfuncs);
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
xlswrite(xlsfilename, errmean, 'Error', 'A2:A29');
xlswrite(xlsfilename, {'St. D.'}, 'Error', 'B1');
xlswrite(xlsfilename, errstd, 'Error', 'B2:B29');
xlswrite(xlsfilename, {'SR'}, 'Error', 'C1');
xlswrite(xlsfilename, succrate, 'Error', 'C2:C29');
xlswrite(xlsfilename, G, 'Convergence', 'A1:U1');
xlswrite(xlsfilename, errmedian, 'Convergence', 'A2:U29');

% Convergence Graph (Example)
figure;
semilogy(G, errmedian(14, :), 'k');
title('Function 14');
xlabel('Generation');
ylabel('Solution Error');

% Dynamic of q value (Example)
figure;
boxplot(q(:, :, 14), G, 'colors', 'k', 'plotstyle','compact');
title('Function 14');
xlabel('Generation');
ylabel('q');

% Dynamic of MF value (Example)
if isfield(allout{1, 1}, 'MF')
	figure;
	boxplot(MF(:, :, 14), G, 'colors', 'k', 'plotstyle','compact');
	title('Function 14');
	xlabel('Generation');
	ylabel('MF');
end

% Dynamic of MCR value (Example)
if isfield(allout{1, 1}, 'MCR')
	figure;
	boxplot(MCR(:, :, 14), G, 'colors', 'k', 'plotstyle','compact');
	title('Function 14');
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