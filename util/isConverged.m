function ret = isConverged(X, TolX)
stdX = std(X, 0, 2);
meanX = mean(abs(X), 2);
ret = all(stdX <= 10 * eps(meanX)) || ...
	all(stdX <= 10 * realmin) || all(stdX <= TolX);
end
