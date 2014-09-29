measureOptions.Dimension = 50;
solvers = {...
	'lshade_sps'};
NP = 18 * measureOptions.Dimension;
Q = [16, 32, 64];
NPmin = {'4'};
F = 0.5;
R = 0;
solverOptions.CR = 0.5;
H = 6;

for i = 1 : numel(NP)
	for j = 1 : numel(Q)
		for k = 1 : numel(NPmin)
			for m = 1 : numel(F)
				for n = 1 : numel(R)
					for p = 1 : numel(H)
						solverOptions.NP = NP(i);
						solverOptions.NPmin = NPmin{k};
						solverOptions.F = F(m);
						solverOptions.R = R(n);
						solverOptions.H = H(p);
						run_complete_cec14(...
							solvers, ...
							Q(j), ...
							solverOptions, ...
							measureOptions);
					end
				end
			end
		end
	end
end
