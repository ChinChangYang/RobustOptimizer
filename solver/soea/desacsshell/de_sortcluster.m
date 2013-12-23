function c = de_sortcluster(f, K, c)
fcluster = zeros(1, K);

for i = 1 : K
	fci = f(c == i);
	fcluster(i) = mean(fci);
end

[~, index] = sort(fcluster);

for i = 1 : length(c)
	c(i) = index(c(i));
end

end