function ret = ranksumtest(X, Y)
[~, N] = size(X);
ret = char(1, N);

for i = 1 : N
	[~, h, stats] = ranksum(X(:, i), Y(:, i), 'method', 'approximate');
	if h == 0
		ret(i) = '=';
	elseif stats.zval < 0
		ret(i) = '-';
	elseif stats.zval > 0
		ret(i) = '+';
	else
		ret(i) = '?';
	end
end
end

