solvers = {...
	'lshade'};
NP = 540;
Q = 64;
NPmin = {'4'};
beta = 1;
R = 0;
solverOptions.F = 0.5;
solverOptions.CR = 0.5;

for i = 1 : numel(NP)
	for j = 1 : numel(Q)
		for k = 1 : numel(NPmin)
			for m = 1 : numel(beta)
				for n = 1 : numel(R)
					solverOptions.NP = NP(i);
					solverOptions.NPmin = NPmin{k};
					solverOptions.R = R(n);
					run_complete_cec14(solvers, Q(j), solverOptions);
				end
			end
		end
	end
end
