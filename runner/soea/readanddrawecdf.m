clear;
%% Original solution errors/function values
eo = zeros(51, 8, 94);

load('SPSDE_CEC11/filenames_o_201405241045.mat');
filenames_tmp = filenames;
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC11/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 4);
	eo(:, i, 1 : 4) = ei;
end

load('SPSDE_CEC14_D10/filenames_201405160848.mat');
filenames_tmp = filenames;
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC14_D10/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 30);
	eo(:, i, 4 + 1 : 4 + 30) = ei;
end

load('SPSDE_CEC14_D30/filenames_201405130927.mat');
filenames_tmp = filenames;
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC14_D30/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 30);
	eo(:, i, 1 + 4 + 30 : 4 + 30 + 30) = ei;
end

load('SPSDE_CEC14_D50/filenames_201405141315.mat');
filenames_tmp = filenames;
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC14_D50/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 30);
	eo(:, i, 1 + 4 + 30 + 30 : 4 + 30 + 30 + 30) = ei;
end

%% Solution errors/function values with the SPS frameworks
esps = zeros(51, 8, 94);

load('SPSDE_CEC11/filenames_sps_201405241045.mat');
filenames_tmp = filenames;
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC11/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 4);
	esps(:, i, 1 : 4) = ei;
end

load('SPSDE_CEC14_D10/filenames_201405231021.mat');
filenames_tmp = filenames(6, :);
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC14_D10/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 30);
	esps(:, i, 4 + 1 : 4 + 30) = ei;
end

load('SPSDE_CEC14_D30/filenames_201405251348.mat');
filenames_tmp = filenames(6, :);
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC14_D30/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 30);
	esps(:, i, 1 + 4 + 30 : 4 + 30 + 30) = ei;
end

load('SPSDE_CEC14_D50/filenames_201405261553.mat');
filenames_tmp = filenames(6, :);
for i = 1 : numel(filenames_tmp)
	filename = sprintf('SPSDE_CEC14_D50/%s', filenames_tmp{i});
	load(filename);
	ei = reshape(allfvals(end, :, :), 51, 1, 30);
	esps(:, i, 1 + 4 + 30 + 30 : 4 + 30 + 30 + 30) = ei;
end
save('ecdfraw.mat', 'eo', 'esps');

%% Compute ECDF
[nruns, nA, nf] = size(eo);
eomin = reshape(min(min(eo)), nf, 1);
eomax = reshape(max(max(eo)), nf, 1);
espsmin = reshape(min(min(esps)), nf, 1);
espsmax = reshape(max(max(esps)), nf, 1);
emin = min(eomin, espsmin);
emax = max(eomax, espsmax);
eonorm = eo;
espsnorm = esps;

for i = 1 : nruns
	for j = 1 : nA
		for k = 1 : nf
			eonorm(i, j, k) = ... 
				(eo(i, j, k) - emin(k) + eps) ./ (emax(k) - emin(k) + eps);
			
			espsnorm(i, j, k) = ... 
				(esps(i, j, k) - emin(k) + eps) ./ (emax(k) - emin(k) + eps);
		end
	end
end

h1 = figure;
hold off;
plot(sort(eonorm(:)), (1 : numel(eonorm)) ./ numel(eonorm), 'b');
hold on;
plot(sort(espsnorm(:)), (1 : numel(espsnorm)) ./ numel(espsnorm), 'r--');
title('ECDF over four real-world problems and 90 artificial functions');
xlabel('NSE');
ylabel('ECDF over all functions');
legend('DEs without the proposed framework', ...
	'DEs with the proposed framework', ...
	'Location', 'SouthEast');
h1Position = get(h1, 'Position');
set(h1, 'Position', [h1Position(1:2), 400 * 1.5, 320 * 1.5]);