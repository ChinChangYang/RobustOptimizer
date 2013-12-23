function run_ranksum_test
% RUN_RANKSUM_TEST Run rank sum test on the data in mat files
addprojectpath;
matfiles = {'err_minmaxcodebest1bin', ...
	'err_minmaxtcde', 'err_minmaxdegl'};

for i = 1 : numel(matfiles)
	basematfile = matfiles{i};
	load(basematfile);
	base_err = err;
	
	for j = 1 : numel(matfiles)
		targetmatfile = matfiles{j};
		load(targetmatfile);
		target_err = err;
		[nFunctions, nMaxFunEvals, ~] = size(target_err);
		results = cell(nFunctions, nMaxFunEvals);
		countEq = 0;
		countPos = 0;
		countNeg = 0;
		
		for k = 1 : nFunctions
			for m = 1 : nMaxFunEvals
				base_err_km = base_err(k, m, :);
				target_err_km = target_err(k, m, :);
				[~, h, stats] = ranksum(base_err_km(:), target_err_km(:));
				
				if h == 0
					results{k, m} = '=';
					countEq = countEq + 1;
				elseif stats.zval >= 0
					results{k, m} = '+';
					countPos = countPos + 1;
				else
					results{k, m} = '-';
					countNeg = countNeg + 1;
				end
				
				fprintf('Base: %s, Target: %s, Func: #%d, MaxFunEvals: #%d, Sign: %s\n', ...
					basematfile, targetmatfile, k, m, results{k, m});
			end
		end
	end
end
end
