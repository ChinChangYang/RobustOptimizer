if matlabpool('size') == 0
	matlabpool('open');
end

clear;
close all;
solvers = {'derand1bin_s'};

for isolver = 1 : numel(solvers)
	startTime = tic;
	date = datestr(now, 'yyyymmddHHMM');
	
	solver = solvers{isolver};
	measureOptions.Dimension = 30;
	measureOptions.Runs = 51;
	measureOptions.FitnessFunctions = {...
		'cec13_f1', 'cec13_f2', 'cec13_f3', 'cec13_f4', ...
		'cec13_f5', 'cec13_f6', 'cec13_f7', 'cec13_f8', 'cec13_f9', ...
		'cec13_f10', 'cec13_f11', 'cec13_f12', 'cec13_f13', ...
		'cec13_f14', 'cec13_f15', 'cec13_f16', 'cec13_f17', ...
		'cec13_f18', 'cec13_f19', 'cec13_f20', 'cec13_f21', ...
		'cec13_f22', 'cec13_f23', 'cec13_f24', 'cec13_f25', ...
		'cec13_f26', 'cec13_f27', 'cec13_f28'};
	measureOptions.MaxFunEvalSet = measureOptions.Dimension * 1e4;
	solverOptions.dimensionFactor = 5;
	solverOptions.NP = solverOptions.dimensionFactor * measureOptions.Dimension;
	solverOptions.F = 0.7;
	solverOptions.CR = 0.5;
	solverOptions.RecordPoint = 20;
	solverOptions.ftarget = 1e-8;
	solverOptions.TolStagnationIteration = inf;
	[allerr, allfes] = err_fes_cec13(solver, measureOptions, solverOptions);
	[nmaxfunevals, nfunctions, nruns] = size(allerr);
	succ_perform = zeros(nmaxfunevals, nfunctions);
	for i = 1 : nmaxfunevals
		for j = 1 : nfunctions
			succ_indices = allerr(i, j, :) <= 1e-8;
			if sum(succ_indices(:)) > 0
				succ_perform(i, j) = ...
					mean(allfes(succ_indices)) / mean(succ_indices(:));
			else
				succ_perform(i, j) = inf;
			end
		end
	end
	
	elapsed_time = toc(startTime);
	if elapsed_time < 60
		fprintf('Elapsed time is %f seconds\n', elapsed_time);
	elseif elapsed_time < 60*60
		fprintf('Elapsed time is %f minutes\n', elapsed_time/60);
	elseif elapsed_time < 60*60*24
		fprintf('Elapsed time is %f hours\n', elapsed_time/60/60);
	else
		fprintf('Elapsed time is %f days\n', elapsed_time/60/60/24);
	end
	
	save(sprintf('succ_perform_cec13_%s.mat', date), ...
		'succ_perform', ...
		'allerr', ...
		'allfes', ...
		'solver', 'measureOptions', 'solverOptions', 'elapsed_time');
end

beep; pause(0.7);
beep; pause(0.31);
beep; pause(0.31);
