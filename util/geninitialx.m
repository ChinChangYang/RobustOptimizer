rng('default');
lb = -100;
ub = 100;
D = [10, 20, 30, 50, 100];
for i = 1 : numel(D)
	Di = D(i);
	NP = 5 * Di;
	X = ...
		repmat(lb, Di, NP) + ...
		repmat(ub - lb, Di, NP) .* lhsdesign(NP, Di, 'iteration', 100)';
	Xname = sprintf('XD%dNP%d', Di, NP);
	eval(sprintf('%s = X;', Xname));
	save('InitialX.mat', Xname, '-append');
end
