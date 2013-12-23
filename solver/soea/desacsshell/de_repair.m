function u = de_repair(u, opts)
%DE_REPAIR Repair of Differential Evolution

if opts.de_handleboundary
	[D, NP] = size(u);
	lb = opts.de_lb;
	ub = opts.de_ub;
	
	if length(lb) < D
		lb = repmat(lb(1), D, 1);
	elseif length(lb) > D
		lb = lb(1 : D);
	end
	
	if length(ub) < D
		ub = repmat(ub(1), D, 1);
	elseif length(ub) > D
		ub = ub(1 : D);
	end
	
	for i = 1 : NP
		for j = 1 : D
			if u(j, i) < lb(j)
				u(j, i) = 2 * lb(j) - u(j, i);
			elseif u(j, i) > ub(j)
				u(j, i) = 2 * ub(j) - u(j, i);
			end
		end
	end
end
end

