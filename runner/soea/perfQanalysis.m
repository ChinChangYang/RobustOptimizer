clear;
load('filenames_201404182025.mat');
[nQ, nA] = size(filenames);
filenames_tmp = filenames;
load(filenames_tmp{1, 1});
[~, ~, nF] = size(allfvals);
ME = zeros(nF, nA, nQ);

for i = 1 : nQ
	for j = 1 : nA
		load(filenames_tmp{i, j});
		ME(:, j, i) = mean(allfvals(end, :, :), 2);
	end
end

MEmin = min(min(ME, [], 2), [], 3);
MEmax = max(max(ME, [], 2), [], 3);
NME = (ME - repmat(MEmin, [1, nA, nQ])) ./ ...
	repmat(MEmax - MEmin + realmin, [1, nA, nQ]);
MNME = sum(sum(NME), 2) / nF / nA;
MNME = MNME(:);

close all;
semilogx(Q, MNME, 'kx-');
title('MNME over CEC 2014 functions');
xlabel('Q');
ylabel('MNME');
