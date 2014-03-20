if matlabpool('size') == 0
    matlabpool('open');
end

clear;
close all;
startTime = tic;
date = datestr(now, 'yyyymmddHHMM');

solver = 'jadebin';
measureOptions.Dimension = 10;
measureOptions.Runs = 51;
measureOptions.FitnessFunctions = {...
	'cec14_f1', 'cec14_f2', 'cec14_f3', 'cec14_f4', ...
	'cec14_f5', 'cec14_f6', 'cec14_f7', 'cec14_f8', 'cec14_f9', ...
	'cec14_f10', 'cec14_f11', 'cec14_f12', 'cec14_f13', ...
	'cec14_f14', 'cec14_f15', 'cec14_f16', 'cec14_f17', ...
	'cec14_f18', 'cec14_f19', 'cec14_f20', 'cec14_f21', ...
	'cec14_f22', 'cec14_f23', 'cec14_f24', 'cec14_f25', ...
	'cec14_f26', 'cec14_f27', 'cec14_f28', 'cec14_f29', ...
	'cec14_f30'};
measureOptions.MaxFunEvalSet = measureOptions.Dimension * 1e4;
solverOptions.RecordPoint = 0;
results = errors_cec2014(solver, measureOptions, solverOptions);

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

save(sprintf('errors_cec2014_%s.mat', date), ...
	'results', 'solver', 'measureOptions', 'solverOptions', 'elapsed_time');