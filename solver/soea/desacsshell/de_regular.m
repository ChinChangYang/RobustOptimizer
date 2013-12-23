function [x, f, counteval] = de_regular(x, K, c, opts, f, counteval, fitfun, varargin)
for i = 1 : K
	[Dim, xciNP] = size(x(:, c == i));
	if xciNP >= opts.de_minlocalpopsize
		if cond(x(:, c == i)) > 1e12
			[V, D] = eig(cov(x(:, c == i)'));
			D = D + 1e-12 * eye(Dim);
			indexes = find(c == i);
			for j = indexes'				
				x(:, j) = real(x(:, j) + 1e-9 * x(:, j) .* (V * sqrt(D) * randn(Dim, 1)));
				f(j) = feval(fitfun, x(:, j), varargin{:});
				counteval = counteval + 1;				
			end
		end
	end
end
end