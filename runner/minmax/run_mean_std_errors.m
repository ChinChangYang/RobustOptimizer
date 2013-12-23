function run_mean_std_errors
%RUN_MEAN_STD_ERRORS Compute the mean and standard deviation of the
%solution errors in mat files
matfiles = {'err_minmaxtcderand1bin', 'err_minmaxtcdebest1bin', ...
	'err_minmaxtcjadebin', 'err_minmaxdegl'};

for i = 1 : numel(matfiles)
	matfile = matfiles{i};
	load(matfile);	
	[nFunctions, nMaxFunEvals, ~] = size(err); %#ok<NODEF>
	fprintf('Data: %s\n', matfile);
	for j = 1 : nFunctions
		for k = 1 : nMaxFunEvals
			err_jk = err(j, k, :);
			fprintf('%.2E\t%.2E\n', ...
				mean(err_jk(:)), std(err_jk(:)));
		end
	end
end

for i = 1 : numel(matfiles)
	matfile = matfiles{i};
	load(matfile);	
	[nFunctions, nMaxFunEvals, ~] = size(err); 
	for j = 1 : nFunctions
		for k = 1 : nMaxFunEvals
			err_jk = err(j, k, :);
			fprintf('Data: %s, Func: #%d, MaxFunEvals: #%d\n', ...
				matfile, j, k);
			fprintf('Mean: %.2E, St. D.: %.2E\n', ...
				mean(err_jk(:)), std(err_jk(:)));
		end
	end
end
end

