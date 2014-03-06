% function analyze_np
if matlabpool('size') == 0
    matlabpool('open');
end

clear;
close all;
startTime = tic;
date = datestr(now, 'yyyymmddHHMM');
RUN = 11;
NP = [4,6,9,14,21,32,48,72,108,162,243,365,548,822,1233];
% NP = [4, 1233];
% D = [10, 30, 50];
D = 30;
fitfun = {'cec13_f1', 'cec13_f2', 'cec13_f3', 'cec13_f4', ...
		'cec13_f5', 'cec13_f6', 'cec13_f7', 'cec13_f8', 'cec13_f9', ...
		'cec13_f10', 'cec13_f11', 'cec13_f12', 'cec13_f13', ...
		'cec13_f14', 'cec13_f15', 'cec13_f16', 'cec13_f17', ...
		'cec13_f18', 'cec13_f19', 'cec13_f20', 'cec13_f21', ...
		'cec13_f22', 'cec13_f23', 'cec13_f24', 'cec13_f25', ...
		'cec13_f26', 'cec13_f27', 'cec13_f28'};
results = zeros(RUN, numel(fitfun), numel(NP), numel(D));
solver = 'derand1bin';
solverOptions.nonlcon = [];
solverOptions.F = 0.7;
solverOptions.CR = 0.5;
solverOptions.TolX = 0;
solverOptions.TolFun = 0;
solverOptions.TolStagnationIteration = 60;
solverOptions.ftarget = -Inf;
solverOptions.Display = 'off';
solverOptions.RecordPoint = 1000;

for Di = 1 : numel(D)
	lb = -100 * ones(D(Di), 1);
	ub = 100 * ones(D(Di), 1);
	maxfunevals = D(Di) * 1e4;
	for NPi = 1 : numel(NP)
		solverOptions.dimensionFactor = NP(NPi)/D(Di);
		for Fi = 1 : numel(fitfun)
			fitfuni = fitfun{Fi};
			parfor RUNi = 1 : RUN
				[~, fmin, ~] = ...
					feval(solver, fitfuni, lb, ub, maxfunevals, solverOptions);
				if fmin <= 1e-8
					fmin = 1e-8;
				end
				results(RUNi, Fi, NPi, Di) = fmin;
			end
		end		
		
		save(sprintf('analyze_np_%s.mat', date), ...
			'results');
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
			fprintf('Min solution error = %.4E\n', min_results(NPi, Di));
			fprintf('25%% quantile of solution error = %.4E\n', quantile_results_25(NPi, Di));
			fprintf('Median solution error = %.4E\n', median_results(NPi, Di));
			fprintf('75%% quantile of solution error = %.4E\n', quantile_results_75(NPi, Di));
			fprintf('Max solution error = %.4E\n', max_results(NPi, Di));
		end
	end
end

save(sprintf('analyze_np_%s.mat', date), ...
	'results');

normalized_results = results;
for RUNi = 1 : RUN
	for Fi = 1 : numel(fitfun)
		results_Fi = results(:, Fi, :);		
		for NPi = 1 : numel(NP)
			normalized_results(RUNi, Fi, NPi) = ...
				(1 - 1e-8) * ...
				(results(RUNi, Fi, NPi) - min(results_Fi(:)) + 1e-8) ./ ...
				(max(results_Fi(:)) - min(results_Fi(:)) + 1e-8) + ...
				1e-8;
		end
	end
end

compact_nr = reshape(normalized_results, RUN * numel(fitfun), numel(NP));
min_nr = reshape(min(compact_nr), numel(NP), 1);
quantile25_nr = reshape(quantile(compact_nr, 0.25), numel(NP), 1);
median_nr = reshape(median(compact_nr), numel(NP), 1);
quantile75_nr = reshape(quantile(compact_nr, 0.75), numel(NP), 1);
max_nr = reshape(max(compact_nr), numel(NP), 1);

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
% end