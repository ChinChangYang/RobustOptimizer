if matlabpool('size') == 0
    matlabpool('open');
end

clear;
close all;
startTime = tic;
date = datestr(now, 'yyyymmddHHMM');

solver = 'mshadeeig';
measureOptions.Dimension = 30;
measureOptions.Runs = 51;
measureOptions.FitnessFunctions = {...
	'cec13_f1', 'cec13_f2', 'cec13_f3', 'cec13_f4', ...
	'cec13_f5'};
measureOptions.MaxFunEvalSet = measureOptions.Dimension * 1e4;
solverOptions.RecordPoint = 0;
solverOptoins.cc = 0.05;
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

beep; pause(0.1);
beep; pause(0.1);
beep; pause(0.1);
