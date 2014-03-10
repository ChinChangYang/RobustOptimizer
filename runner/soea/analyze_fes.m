if matlabpool('size') == 0
    matlabpool('open');
end

clear;
close all;
startTime = tic;
date = datestr(now, 'yyyymmddHHMM');
RUN = 4;
NP = [4,6,9,14,21,32,48,72,108,162,243,365,548,822,1233,1850,2775,4163];
% NP = [4, 4163];
% D = [10, 30, 50];
D = 50;
fitfun = {'cec13_f1'};
results = zeros(RUN, numel(fitfun), numel(NP), numel(D));
solver = 'derand1bin';
solverOptions.nonlcon = [];
solverOptions.F = 0.7;
solverOptions.CR = 0.5;
solverOptions.TolX = 1;
solverOptions.TolFun = 0;
solverOptions.TolStagnationIteration = 60;
solverOptions.ftarget = -Inf;
solverOptions.Display = 'off';
solverOptions.RecordPoint = 0;

for Di = 1 : numel(D)
	lb = -100 * ones(D(Di), 1);
	ub = 100 * ones(D(Di), 1);
	maxfunevals = D(Di) * 1e5;
	for NPi = 1 : numel(NP)
		solverOptions.dimensionFactor = NP(NPi)/D(Di);
		for Fi = 1 : numel(fitfun)
			fitfuni = fitfun{Fi};
			fprintf('NP: %d; fitfun: %s\n', NP(NPi), fitfuni);
			parfor RUNi = 1 : RUN
				[~, ~, out] = ...
					feval(solver, fitfuni, lb, ub, maxfunevals, solverOptions);

				results(RUNi, Fi, NPi, Di) = out.fes;
			end
		end		
		
		save(sprintf('analyze_fes_%s.mat', date), ...
			'results', 'D', 'NP', 'RUN', 'fitfun');
	end
end

min_results = reshape(min(results), numel(fitfun), numel(NP), numel(D));
quantile_results_25 = reshape(quantile(results, 0.25), numel(fitfun), numel(NP), numel(D));
median_results = reshape(median(results), numel(fitfun), numel(NP), numel(D));
quantile_results_75 = reshape(quantile(results, 0.75), numel(fitfun), numel(NP), numel(D));
max_results = reshape(max(results), numel(fitfun), numel(NP), numel(D));

for Di = 1 : numel(D)
	for Fi = 1 : numel(fitfun)
		for NPi = 1 : numel(NP)
			fprintf('D = %d; NP = %d; fitfun = %s\n', D(Di), NP(NPi), fitfun{Fi});
			fprintf('Min FEs = %.4E\n', min_results(Fi, NPi, Di));
			fprintf('25%% quantile of FEs = %.4E\n', quantile_results_25(Fi, NPi, Di));
			fprintf('Median FEs = %.4E\n', median_results(Fi, NPi, Di));
			fprintf('75%% quantile of FEs = %.4E\n', quantile_results_75(Fi, NPi, Di));
			fprintf('Max FEs = %.4E\n', max_results(Fi, NPi, Di));
		end
	end
end

min_results = min_results';
quantile_results_25 = quantile_results_25';
median_results = median_results';
quantile_results_75 = quantile_results_75';
max_results = max_results';

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

save(sprintf('analyze_fes_%s.mat', date), ...
	'results', 'D', 'NP', 'RUN', 'fitfun', 'elapsed_time');
% end