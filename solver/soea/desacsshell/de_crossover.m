function u = de_crossover(x, v, c, K, opts)
%DE_CROSSOVER Crossover of Differential Evolution

CR = opts.de_cr;
[D, xNP] = size(x);
u = v;
MAXK = opts.de_maxmutatecluster;

if K > MAXK
	K = ceil(K / 2);
end

for i = 1 : K
	xci = x(:, c == i);
	[~, xci_NP] = size(x(:, c == i));
	
	if xci_NP > 1
		[B, ~] = eig(cov(xci'));
		
		for j = 1 : ceil(xNP / K)
			jrand = floor(1 + D * rand);
			r1 = floor(1 + xci_NP * rand);
			
			if rand < 0.1
				Bjrand = B(:, jrand);
				index = (i - 1) * ceil(xNP / K) + j;
				L = dot(v(:, index) - xci(:, r1), Bjrand);
				if L < 1e-12 && L > -1e-12
					L = L + 1e-12 * (4 + randn);
				end
					
				u(:, index) = ...
					xci(:, r1) + L * Bjrand;
			end
		end
	end
end

end

