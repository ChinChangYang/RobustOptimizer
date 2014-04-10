fileQ1 = {'cec13D30_dcmaea_s_Q10_201404101119.mat', ...
	'cec13D30_dcmaea_s_Q20_201404101127', ...
	'cec13D30_dcmaea_s_Q30_201404101134'};
fileQ2 = {'cec13D30_dcmaea_s_Q40_201404101142', ...
	'cec13D30_dcmaea_s_Q30_201404101134', ...
	'cec13D30_dcmaea_s_Q50_201404101149'};
fileQ3 = {'cec13D30_dcmaea_s_Q40_201404101142', ...
	'cec13D30_dcmaea_s_Q50_201404101149', ...
	'cec13D30_dcmaea_s_Q60_201404101157'};
fileQall = [fileQ1; fileQ2; fileQ3];

[nQ, nA] = size(fileQall);
load(fileQall{1, 1});
[~, ~, nF] = size(allfvals);
MNME = zeros(nQ, 1);

for i = 1 : nQ
	ME = zeros(nF, nA);
	for j = 1 : nA
		load(fileQall{i, j});
		ME(:, j) = mean(allfvals(end, :, :), 2);
	end
	MEmin = min(ME, [], 2);
	MEmax = max(ME, [], 2);
	NME = (ME - repmat(MEmin, 1, nA) + eps) ./ repmat(MEmax - MEmin + eps, 1, nA);
	MNME(i) = mean(NME(:));
end
plot(MNME, 'kx-');
title('MNME over CEC 2013 functions');
xlabel('Q');
ylabel('MNME');
