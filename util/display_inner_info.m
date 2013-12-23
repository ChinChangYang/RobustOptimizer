function display_inner_info(innerState)
%DISPLAY_INNER_INFO Display information of inner states
n_innerState = numel(innerState);

noNaN_index = -1;
for i = 1 : n_innerState
	if ~isempty(innerState{i})
		noNaN_index = i;
		break;
	end
end

[D, ~] = size(innerState{noNaN_index}.X);
stdX = nan(D, n_innerState);
fstd = nan(1, n_innerState);
N_noNaN = 0;

for i = 1 : n_innerState
	if ~isempty(innerState{i})
		stdX(:, i) = std(innerState{i}.X, 0, 2);
		fstd(i) = std(innerState{i}.f);
		N_noNaN = N_noNaN + 1;
	end
end

fstd(isnan(fstd)) = [];
stdX_noNaN = zeros(D, N_noNaN);
i_stdX_noNaN = 1;

for i = 1 : n_innerState
	if all(~isnan(stdX(:, i)))
		stdX_noNaN(:, i_stdX_noNaN) = stdX(:, i);
		i_stdX_noNaN = i_stdX_noNaN + 1;
	end
end

stdX = stdX_noNaN;

minStdX = min(stdX, [], 2);
maxStdX = max(stdX, [], 2);
meanStdX = mean(stdX(:));
medianStdX = median(stdX(:));

[maxmaxStdX, maxStdXIdx] = max(maxStdX);
[minminStdX, minStdXIdx] = min(minStdX);
fprintf(['InnerState: ', ...
	'Fstd: %0.2E,\t', ...
	'Min Xstd(%d): %0.2E,\t', ...
	'Max Xstd(%d): %0.2E,\t', ...
	'Mean Xstd: %0.2E,\t', ...
	'Median Xstd: %0.2E\n', ...
	'--\n'], ...
	mean(fstd), ...
	minStdXIdx, minminStdX, ...
	maxStdXIdx, maxmaxStdX, ...
	meanStdX, medianStdX);

end

