function v = de_mutate(x, K, NP, F, c, opts)
%DE_MUTATE Mutation of Differential Evolution

minlocalpopsize = opts.de_minlocalpopsize;
PF = opts.de_pf;
MAXK = opts.de_maxmutatecluster;
v = x;
[~, xNP] = size(x);

if K > MAXK
	K = ceil(K / 2);
end

for i = 1 : K
	xciNP = sum(c == i);
	if xciNP >= minlocalpopsize
		xci = x(:, c == i);
		for j = 1 : ceil(NP / K)
			r1 = floor(1 + xciNP * rand);
			r2 = floor(1 + xciNP * rand);
			r3 = r2;
			while r2 == r3
				r3 = floor(1 + xciNP * rand);
			end
			
			v(:, (i - 1) * ceil(NP / K) + j) = ...
				xci(:, r1) + F * (1 + PF * randn) * (xci(:, r2) - xci(:, r3));
		end
	else
		for j = 1 : ceil(NP / K)
			r1 = floor(1 + xNP * rand);
			r2 = floor(1 + xNP * rand);
			r3 = r2;
			while r2 == r3
				r3 = floor(1 + xNP * rand);
			end
			
			v(:, (i - 1) * ceil(NP / K) + j) = ...
				x(:, r1) + F * (1 + PF * randn) * (x(:, r2) - x(:, r3));
		end
	end
end
end

