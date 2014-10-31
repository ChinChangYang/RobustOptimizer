function run_complete_cec14_all
close all;
measureOptions.Dimension = 10;
measureOptions.Runs = 28;
measureOptions.MaxFunEvals = measureOptions.Dimension * 1e4;

solver = 'umoeas_a';
NP1 = 15 .* measureOptions.Dimension;
NP2 = 5 .* measureOptions.Dimension;
Q = 64;
NPmin = {'4'};
F = 0.5;
H = [12, 24, 48];
cw = [0.025, 0.05, 0.1];
CRmax = 0.3;
CRmin = 0.05;
CS = [25, 50, 100];
solverOptions.CR = 0.5;

if matlabpool('size') == 0
	matlabpool('open');
end

filenames = cell(...
	numel(NP1) * numel(NP2) * numel(Q) * numel(NPmin) * numel(F) * ...
	numel(H) * numel(CRmax) * numel(cw) * numel(CRmin), 1);

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
	
	metafiledate = datestr(now, 'yyyymmddHHMM');
	metafilename = sprintf('filenames_%s.mat', metafiledate);
	save(metafilename, 'filenames');
end