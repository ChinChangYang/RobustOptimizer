function simulate_distribution_issue11
T = 1e4;
w = [1.73858, 7, 3.00003, 1.72803, 3.26699];
s = [0, 1, 1.81818, 4, 5];
mu = [0.9091, 1.00001, 1, 1, 2.35795];
w = w / sum(w);
r = [];
for i = 1 : numel(w)
	wT = floor(w(i) * T);
	r = [r, s(i) + exprnd(mu(i), 1, wT)]; %#ok<AGROW>
end

hist(r, 100);
fprintf('Sample Mean: %0.4E\n', mean(r));
fprintf('Expected Mean: %0.4E\n', sum(w .* (s + mu)));
end