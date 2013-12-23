function ret = computeInnerFstd( innerState )
%COMPUTEINNERFSTD Compute standard deviation of function values in inner
%states

n_innerState = numel(innerState);
fstd = nan(1, n_innerState);

for i = 1 : n_innerState
	if ~isempty(innerState{i})
		fstd(i) = std(innerState{i}.f);
	end
end

fstd(isnan(fstd)) = [];
ret = mean(fstd);
end

