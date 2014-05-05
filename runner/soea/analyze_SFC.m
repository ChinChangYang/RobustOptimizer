function analyze_SFC
if matlabpool('size') == 0
	matlabpool('open');
end

startTime = tic;

load('InitialX.mat');
filename = 'analyze_SFC.xlsx';
solver = 'debest1bin_e';
measureOptions.Dimension = 30;
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
lb = -100 * ones(measureOptions.Dimension, 1);
ub = 100 * ones(measureOptions.Dimension, 1);
solverOptions.NP = 5 * measureOptions.Dimension;
solverOptions.ftarget = 1e-8;
solverOptions.Display = 'off';
solverOptions.RecordPoint = 21;
solverOptions.initial.X = eval(sprintf('XD%dNP%d', ...
	measureOptions.Dimension, ...
	solverOptions.NP));

title = {'function', 'mean(out.mSFC)', 'mean(out.mFC)'};
range = 'A1:C1';
xlswrite(filename, title, solver, range);
for iFun = 1 : numel(measureOptions.FitnessFunctions)
	rng('default');
	fitfun = measureOptions.FitnessFunctions{iFun};
	[~, ~, out] = ...
		feval(solver, fitfun, lb, ub, measureOptions.MaxFunEvalSet, solverOptions);
	
	range = sprintf('A%d', iFun + 1);
	xlswrite(filename, iFun, solver, range);
	range = sprintf('B%d:C%d', iFun + 1, iFun + 1);
	xlswrite(filename, [mean(out.mSFC), mean(out.mFC)], solver, range);
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
end

