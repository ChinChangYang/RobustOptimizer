solvers = {...
	'shade_sps_eig'};
NP = [60, 120, 180];
Q = [32, 64, 128];
R = [0.3, 0.5, 0.7];

for i = 1 : numel(NP)
	for j = 1 : numel(Q)
		for k = 1 : numel(R)
			run_complete_cec14(solvers, NP(i), Q(j), R(k));
		end
	end
end
