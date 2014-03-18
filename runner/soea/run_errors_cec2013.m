if matlabpool('size') == 0
    matlabpool('open');
end

clear;
close all;
startTime = tic;
date = datestr(now, 'yyyymmddHHMM');

solver = 'jadebin';
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
solverOptions.RecordPoint = 0;
results = errors_cec2013(solver, measureOptions, solverOptions);

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

save(sprintf('errors_cec2013_%s.mat', date), ...
	'results', 'solver', 'measureOptions', 'solverOptions', 'elapsed_time');