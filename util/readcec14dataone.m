function readcec14dataone(metafilename)
load(metafilename);
filenames_o = filenames;
nA = numel(filenames_o);
load(filenames_o{1});
[~, nf] = size(allout);

D					= measureOptions.Dimension;
tablefilename		= sprintf('CEC14_D%d_TABLE.xlsx', D);
convfilename		= sprintf('CEC14_D%d_CONV.xlsx', D);
contfilename		= sprintf('CEC14_D%d_CONT.xlsx', D);

meansuccrates			= zeros(1, nA);
errmeanall				= zeros(nf, nA);

for i = 1 : nA
	[solver, errmean, errstd, succrate, compcomplex, errmedian, ...
		~, distancemedian, fes, G, ~] = ...
		readonedata(filenames_o{i}); %#ok<ASGLU>
	
	solver = sprintf('%d.%s', i, solver);
	xlswrite(tablefilename, {solver}, solver, 'A1');
	xlswrite(tablefilename, {'Mean'}, solver, 'A2');
	xlswrite(tablefilename, errmean, solver, 'A3:A32');
	xlswrite(tablefilename, {'St. D.'}, solver, 'B2');
	xlswrite(tablefilename, errstd, solver, 'B3:B32');
	xlswrite(tablefilename, {'SR'}, solver, 'C2');
	xlswrite(tablefilename, succrate, solver, 'C3:C32');
	xlswrite(convfilename, G, solver, 'A1:U1');
	xlswrite(convfilename, errmedian, solver, 'A2:U31');
	xlswrite(contfilename, G, solver, 'A1:U1');
	xlswrite(contfilename, distancemedian, solver, 'A2:U31');
	
	meansuccrates(i) = mean(succrate);
	errmeanall(:, i) = errmean;
	fprintf('%d -- Mean Succ. Rate: %.2f%%\n', i, 100 * meansuccrates(i));
end

% Normalize Mean Error
normerrmean = zeros(nf, nA);
for i = 1 : nf
	for j = 1 : nA
		normerrmean(i, j) = ...
			(errmeanall(i, j) - min(errmeanall(i, :)) + eps) ...
			./ (max(errmeanall(i, :)) - min(errmeanall(i, :)) + eps);
	end
end
meanerrmean = mean(normerrmean);

fprintf('-- Sorting --\n');
[meanerrmean, index] = sort(meanerrmean);
meansuccrates = meansuccrates(index);

for i = 1 : nA
	fprintf('%d -- Mean Succ. Rate: %.2f%% (%.2f)\n', ...
		index(i), ...
		100 * meansuccrates(i), ...
		meanerrmean(i));
end

fprintf('%s: OK!\n', tablefilename);
end

function [solver, errmean, errstd, succrate, compcomplex, errmedian, ...
	qmediansum, distancemedian, fes, G, qmedianmean] = readonedata(filename) %#ok<STOUT>

load(filename);

% Generate Measurements
allfvals(allfvals <= 1e-8) = 0; %#ok<NODEF,*SAGROW>
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

if ~isfield(allout{1, 1}, 'FC') %#ok<USENS>
	qmediansum = [];
	qmedianmean = [];
else
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
	% qmedianmax	= max(qmedian, [], 2);
	% qmedianmax	= reshape(qmedianmax, nprogress, nfuncs)';
	qmedian		= reshape(qmedian, NP, nprogress, nfuncs);
	qmediansum	= sum(qmedian > solverOptions.Q);
	qmediansum	= reshape(qmediansum, nprogress, nfuncs)';
end
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