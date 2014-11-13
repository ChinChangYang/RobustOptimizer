function run_complete_cec14_all
close all;

D = [10, 30, 50];

for iD = 1 : 3
	measureOptions.Dimension = D(iD);
	measureOptions.Runs = 28;
	measureOptions.MaxFunEvals = measureOptions.Dimension * 1e4;
	
	% solver = 'umoeas_b';
	% solver = 'lshade_sps';
	% solver = 'lshade_sps_eig_k';
	% solver = 'umoeas_c';
	solver = 'moeas_a';
	NP1 = 19 .* measureOptions.Dimension;
	NP2 = 5 .* measureOptions.Dimension;
	Q = 64;
	NPmin = {'4'};
	F = 0.5;
	H = 6;
	cw = 0.05;
	CRmax = 0.3;
	CRmin = 0.05;
	CS = 50;
	MixStageFactor = 0.1;
	solverOptions.CR = 0.5;
	solverOptions.NP = 19 .* measureOptions.Dimension;
	
	if matlabpool('size') == 0
		matlabpool('open');
	end
	
	filenames = cell(...
		numel(NP1) * numel(NP2) * numel(Q) * numel(NPmin) * numel(F) * ...
		numel(H) * numel(CRmax) * numel(cw) * numel(CRmin) * numel(CS) * ...
		numel(MixStageFactor), 1);
	
	counter = 1;
	for i = 1 : numel(NP1)
		for j = 1 : numel(Q)
			for k = 1 : numel(NPmin)
				for m = 1 : numel(F)
					for p = 1 : numel(H)
						for q = 1 : numel(CRmax)
							for r = 1 : numel(cw)
								for s = 1 : numel(CRmin)
									for t = 1 : numel(NP2)
										for u = 1 : numel(CS)
											for v = 1 : numel(MixStageFactor)
												solverOptions.NP1 = NP1(i);
												solverOptions.NP2 = NP2(t);
												solverOptions.Q = Q(j);
												solverOptions.NPmin = NPmin{k};
												solverOptions.F = F(m);
												solverOptions.H = H(p);
												solverOptions.CRmax = CRmax(q);
												solverOptions.cw = cw(r);
												solverOptions.CRmin = CRmin(s);
												solverOptions.CS = CS(u);
												solverOptions.MixStageFactor = MixStageFactor(v);
												
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
				end
			end
		end
		
		metafiledate = datestr(now, 'yyyymmddHHMM');
		metafilename = sprintf('filenames_%s.mat', metafiledate);
		save(metafilename, 'filenames');
	end
end