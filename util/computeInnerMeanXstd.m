function ret = computeInnerMeanXstd( innerState )
%COMPUTEINNERMEANXSTD Compute mean of standard deviation of solutions X in
%inner states
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
N_noNaN = 0;

for i = 1 : n_innerState
	if ~isempty(innerState{i})
		stdX(:, i) = std(innerState{i}.X, [], 2);
		N_noNaN = N_noNaN + 1;
	end
end

stdX_noNaN = zeros(D, N_noNaN);
i_stdX_noNaN = 1;

for i = 1 : n_innerState
	if all(~isnan(stdX(:, i)))
		stdX_noNaN(:, i_stdX_noNaN) = stdX(:, i);
		i_stdX_noNaN = i_stdX_noNaN + 1;
	end
end

stdX = stdX_noNaN;
ret = mean(stdX(:));
end

