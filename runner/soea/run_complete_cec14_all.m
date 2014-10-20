function run_complete_cec14_all
measureOptions.Dimension = 30;
measureOptions.Runs = 28;
measureOptions.MaxFunEvals = measureOptions.Dimension * 1e4;

solver = 'lshade_sps_eig_g';
NP = [8, 11, 14] .* measureOptions.Dimension;
Q = 64;
NPmin = {'4'};
F = 0.5;
solverOptions.CR = 0.5;
H = 6;
CRmax = [0.2, 0.3, 0.4];
cw = [0.3, 0.4, 0.5];

if matlabpool('size') == 0
	matlabpool('open');
end

filenames = cell(...
	numel(NP) * numel(Q) * numel(NPmin) * numel(F) * numel(H) * ...
	numel(CRmax) * numel(cw), 1);

counter = 1;
for i = 1 : numel(NP)
	for j = 1 : numel(Q)
		for k = 1 : numel(NPmin)
			for m = 1 : numel(F)
				for p = 1 : numel(H)
					for q = 1 : numel(CRmax)
						for r = 1 : numel(cw)
							solverOptions.NP = NP(i);
							solverOptions.Q = Q(j);
							solverOptions.NPmin = NPmin{k};
							solverOptions.F = F(m);
							solverOptions.H = H(p);
							solverOptions.CRmax = CRmax(q);
							solverOptions.cw = cw(r);
							
							innerdate = datestr(now, 'yyyymmddHHMMSS');
							startTime = tic;
							
							[allout, allfvals, allfes, T0, T1, T2] = ...
								complete_cec14(...
								solver, ...
								measureOptions, ...
								solverOptions); %#ok<NASGU,ASGLU>
							
							elapsedTime = toc(startTime); %#ok<NASGU>
														
							filenames{counter} = sprintf('cec14D%d_%s_%s.mat', ...
								measureOptions.Dimension, solver, innerdate);
							
							save(filenames{counter}, ...
								'allout', ...
								'allfvals', ...
								'allfes', ...
								'T0', 'T1', 'T2', ...
								'solver', ...
								'measureOptions', ...
								'solverOptions', ...
								'elapsedTime');
							
							counter = counter + 1;
						end
					end
				end
			end
		end
	end
end

metafiledate = datestr(now, 'yyyymmddHHMM');
metafilename = sprintf('filenames_%s.mat', metafiledate);
save(metafilename, 'filenames');
end