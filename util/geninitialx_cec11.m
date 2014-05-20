rng('default');
D = [6, 20, 26, 22];
for i = 1 : numel(D)
	Di = D(i);
	NP = 5 * Di;
	
	if i == 1
		lb = -6.4 * ones(Di, 1);
		ub = 6.35 * ones(Di, 1);
	elseif i == 2
		lb = zeros(Di, 1);
		ub = 2 * pi * ones(Di, 1);
	elseif i == 3
		[lb, ub] = getlimit_messenger;
	elseif i == 4
		[lb, ub] = getlimit_cassini2;
	else
		error('BUG: Wrong value of i');
	end
	
	X = ...
		repmat(lb, 1, NP) + ...
		repmat(ub - lb, 1, NP) .* lhsdesign(NP, Di, 'iteration', 100)';
	Xname = sprintf('XF%dNP%d', i, NP);
	eval(sprintf('%s = X;', Xname));
	
	if i == 1
		save('InitialX_CEC11.mat', Xname);
	else
		save('InitialX_CEC11.mat', Xname, '-append');
	end
end
