function out = de_record(out, x, f, c, K, NP, F, counteval, countsuccess, opts)
%DE_RECORD Summary of this function goes here
%   Detailed explanation goes here

out.fes(end + 1) = counteval;
out.successrate(end + 1) = countsuccess / NP;
sumstdx = 0;
conditions = ones(1, K);
mindist = inf;
maxdist = 0;

for i = 1 : K
	xci = x(:, c == i);
	sumstdx = sumstdx + std(xci, [], 2);
	[~, xciNP] = size(xci);
	if xciNP >= opts.de_minlocalpopsize
		conditions(i) = cond(xci);
	end
	
	for j = 2 : xciNP
		dist = norm(xci(:, 1) - xci(:, j));
		if dist < mindist
			mindist = dist;
		elseif dist > maxdist
			maxdist = dist;
		end
	end
end

out.stdx(:, end + 1) = sumstdx / K;
out.de_f(end + 1) = F;
out.k(end + 1) = K;
out.popsize(end + 1) = NP;
out.fmean(end + 1) = mean(f);
out.cond(end + 1) = max(conditions);
out.fmin(end + 1) = min(f);
out.mindist(end + 1) = mindist;
out.maxdist(end + 1) = maxdist;
end

