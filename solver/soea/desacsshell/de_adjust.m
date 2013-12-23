function [K, NP, F] = de_adjust(out, D, K, NP, F, counteval, opts)
	T		= round(opts.reliablegeneration);
	alpha	= opts.clusterscale;
	succmin = opts.minimalsuccessrate;
	succmax = opts.maximalsuccessrate;
	
	if length(out.successrate) > T
		temporalsr = out.successrate(end-T:end);
		temporalk = out.k(end-T:end);
		if length(temporalsr(temporalk == K)) > T
			index = find(temporalk == K, T, 'last');
			if mean(temporalsr(index)) < succmin
				if K < ceil(NP / D / 3)
					K = min(ceil(K * alpha), ceil(NP / D / 3));
				else
					F = F / alpha;
				end
			elseif mean(temporalsr(index)) > succmax
				if F < opts.de_f
					K = max(K - 1, 1);
					F = opts.de_f;
				else
					K = max(floor(K / alpha), 1);
				end
				
				if 2 * counteval >= opts.maxfunevals
					NP = max(floor(0.95 * NP), K * 3 * D);
				end
			else
				F = min(F * alpha, opts.de_f);
			end
		end
	end
end