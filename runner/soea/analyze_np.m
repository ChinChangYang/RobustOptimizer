if matlabpool('size') == 0
    matlabpool('open');
end

clear;
close all;
startTime = tic;
date = datestr(now, 'yyyymmddHHMM');
RUN = 20;
NP = [4	5	6	8	10	12	15	19	24	30	37	47	58	73	91	...
	114	142	178	222	278];
% NP = [4	37 278];
% MaxFEs = [3e4, 3e5, 3e6];
% MaxFEs = [3e4, 3e5];
MaxFEs = 3e5;
% D = [10, 30, 50];
D = 30;
fitfun = {'cec13_f1', 'cec13_f2', 'cec13_f3', 'cec13_f4', ...
		'cec13_f5', 'cec13_f6', 'cec13_f7', 'cec13_f8', 'cec13_f9', ...
		'cec13_f10', 'cec13_f11', 'cec13_f12', 'cec13_f13', ...
		'cec13_f14', 'cec13_f15', 'cec13_f16', 'cec13_f17', ...
		'cec13_f18', 'cec13_f19', 'cec13_f20', 'cec13_f21', ...
		'cec13_f22', 'cec13_f23', 'cec13_f24', 'cec13_f25', ...
		'cec13_f26', 'cec13_f27', 'cec13_f28'};
% fitfun = {'cec13_f1', 'cec13_f19'};
error = zeros(RUN, numel(fitfun), numel(NP), numel(MaxFEs), numel(D));
xstd = zeros(RUN, numel(fitfun), numel(NP), numel(MaxFEs), numel(D));
solver = 'mshadeeig';
solverOptions.nonlcon = [];
% solverOptions.F = 0.7;
% solverOptions.CR = 0.5;
solverOptions.TolX = 0;
solverOptions.TolFun = 0;
solverOptions.TolStagnationIteration = 100;
solverOptions.ftarget = -Inf;
solverOptions.Display = 'off';
solverOptions.RecordPoint = 1;

for Di = 1 : numel(D)
	lb = -100 * ones(D(Di), 1);
	ub = 100 * ones(D(Di), 1);
	for Mi = 1 : numel(MaxFEs)
		MaxFEsi = MaxFEs(Mi);
		for NPi = 1 : numel(NP)
			solverOptions.dimensionFactor = NP(NPi)/D(Di);
			solverOptions.NP = NP(NPi);
			for Fi = 1 : numel(fitfun)
				fitfuni = fitfun{Fi};
				fprintf('D: %d, MaxFEs: %.4E, NP: %d; fitfun: %s\n', ...
					D(Di), MaxFEsi, NP(NPi), fitfuni);
				parfor RUNi = 1 : RUN
					[~, fmin, out] = ...
						feval(solver, fitfuni, lb, ub, MaxFEsi, solverOptions);
					if fmin <= 1e-8
						fmin = 1e-8;
					end
					
					error(RUNi, Fi, NPi, Mi, Di) = fmin;
					xstd(RUNi, Fi, NPi, Mi, Di) = mean(out.xstd);
				end
			end
			
			save(sprintf('analyze_np_%s.mat', date), ...
				'error', ...
				'xstd', ...
				'D', ...
				'MaxFEs', ...
				'NP', ...
				'RUN', ...
				'fitfun', ...
				'solver');
		end
	end
end

mean_error = mean(error);
[~, minidx] = min(mean_error, [], 3);
minidx = reshape(minidx, numel(fitfun), numel(MaxFEs), numel(D));
mean_error_bestnp = zeros(numel(fitfun), numel(MaxFEs), numel(D));

mean_xstd = mean(xstd);
mean_xstd_best = min(min(mean_xstd, [], 3), [], 4);
mean_xstd_best = reshape(mean_xstd_best, numel(fitfun), numel(D));
mean_xstd_overall = mean(reshape(xstd, RUN * numel(fitfun), numel(NP), numel(MaxFEs), numel(D)));
mean_xstd_overall = reshape(mean_xstd_overall, numel(NP), numel(MaxFEs), numel(D));

mean_xstd_bestnp = zeros(numel(fitfun), numel(MaxFEs), numel(D));

for i = 1 : numel(fitfun)
	for j = 1 : numel(MaxFEs)
		for k = 1 : numel(D)
			mean_error_bestnp(i, j, k) = ...
				mean_error(1, i, minidx(i, j, k), j, k);
			
			mean_xstd_bestnp(i, j, k) = ...
				mean_xstd(1, i, minidx(i, j, k), j, k);
		end
	end
end

norm_errors = error;
for Di = 1 : numel(D)
	for Mi = 1 : numel(MaxFEs)
		for Fi = 1 : numel(fitfun)
			errors_Fi = error(:, Fi, :, Mi, Di);
			for NPi = 1 : numel(NP)
				for RUNi = 1 : RUN
					norm_errors(RUNi, Fi, NPi, Mi, Di) = ...
						(1 - 1e-8) * ...
						(error(RUNi, Fi, NPi, Mi, Di) - min(errors_Fi(:)) + 1e-8) ./ ...
						(max(errors_Fi(:)) - min(errors_Fi(:)) + 1e-8) + ...
						1e-8;
				end
			end
		end
	end
end

mean_norm_errors = mean(reshape(norm_errors, RUN * numel(fitfun), numel(NP), numel(MaxFEs), numel(D)));

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

save(sprintf('analyze_np_%s.mat', date), ...
	'error', ...
	'xstd', ...
	'D', ...
	'MaxFEs', ...
	'NP', ...
	'RUN', ...
	'fitfun', ...
	'solver', ...
	'elapsed_time');

beep; pause(0.7);
beep; pause(0.31);
beep; pause(0.31);