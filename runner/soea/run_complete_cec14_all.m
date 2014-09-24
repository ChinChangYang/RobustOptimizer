solvers = {...
	'shade_sps_eig_d'};
NP = [200, 300, 400];
Q = 64;
NPmin = {'30', '60', '90'};
beta = [0.5, 0.75, 1];
R = 0.3;

for i = 1 : numel(NP)
	for j = 1 : numel(Q)
		for k = 1 : numel(NPmin)
			for m = 1 : numel(beta)
				for n = 1 : numel(R)
					solverOptions.NP = NP(i);
					solverOptions.NPmin = NPmin{k};
					solverOptions.beta = beta(m);
					solverOptions.H = eval(NPmin{k});
					solverOptions.R = R(n);
					run_complete_cec14(solvers, Q(j), solverOptions);
				end
			end
		end
	end
end
