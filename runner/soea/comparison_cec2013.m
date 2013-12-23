function allErrs = comparison_cec2013
% COMPARISON_CEC2013 Compare solvers in terms of number of win cases and
% ECDF. 
%
% Results:
% Number of win cases (debest1eig): 13
% Number of win cases (derand1eig): 20
% Number of win cases (jadeeig): 26
% Number of win cases (sadeeig): 30
% Func 1, best solver: sadeeig , fval: 1.28E+02
% Func 2, best solver: debest1eig , fval: 4.65E+01
% Func 3, best solver: jadeeig , fval: 2.88E+01
% Func 4, best solver: jadeeig , fval: 2.03E+01
% Func 5, best solver: derand1eig , fval: 1.59E+01
% Func 6, best solver: debest1eig , fval: 9.87E+00
% Func 7, best solver: debest1eig , fval: 5.07E+00
% Func 8, best solver: sadeeig , fval: 2.11E+00
% Func 9, best solver: jadeeig , fval: 1.28E+00
% Func 10, best solver: jadeeig , fval: 3.85E-01
% Func 11, best solver: sadeeig , fval: 2.66E-01
% Func 12, best solver: derand1eig jadeeig , fval: 2.32E-01
% Func 13, best solver: sadeeig , fval: 9.95E-02
% Func 14, best solver: jadeeig , fval: 9.95E-02
% Func 15, best solver: sadeeig , fval: 1.99E-01
% Elapsed time is 99.473145 seconds.
% 
% The above result was gathered with Intel(R) Core(TM) i3-3220 @ 3.30GHz.

comparisonStart = tic;
addprojectpath;

% Determine solvers to be compared
solvers = {'debest1eig', 'derand1eig', 'jadeeig', 'sadeeig'};
nSolvers = numel(solvers);
errs = cell(1, nSolvers);

% Compute errors
D = 2;
measureOptions.Dimension = D;
measureOptions.Runs = 30;
measureOptions.MaxFunEvalSet = round(2.^(linspace(1, 14, 15)));
measureOptions.FitnessFunctions = {'cec13_f12'};
for i = 1 : nSolvers
	solver = solvers(i);
	errs{i} = mean(errors_cec2013(solver{:}, measureOptions), 3);
end

savefilename = sprintf('comparison_cec2013_D%d.mat', D);
save(savefilename, 'errs');

% Compute number of win cases
nWins = zeros(1, nSolvers);
for i = 1 : nSolvers
	for j = 1 : nSolvers
		if i ~= j
			nWins(i) = nWins(i) + sum(errs{i} < errs{j});
		end
	end
end

% Display numbers of win cases for all solvers
for i = 1 : nSolvers
	solver = solvers(i);
	fprintf('Number of win cases (%s): %d\n', solver{:}, nWins(i));
end

% Draw ECDF of the normalized mean errors
hold off;
sAllErrs = '[errs{1}''';
for i = 2 : nSolvers
	sErrs_i = sprintf('errs{%d}''', i);
	sAllErrs = [sAllErrs '; ' sErrs_i]; %#ok<AGROW>
end
sAllErrs = [sAllErrs ']'];
allErrs = eval(sAllErrs);
drawecdf(allErrs);
legend(solvers);

[minErrs, minErrIndexes] = min(allErrs);

for i = 1 : numel(minErrIndexes)
	bestSolverIndexes = find(allErrs(:, i) == minErrs(i));
	bestSolvers = '';
	for j = 1 : numel(bestSolverIndexes)
		jBestSolver = solvers(bestSolverIndexes(j));
		bestSolvers = [bestSolvers, jBestSolver{:}, ' ']; %#ok<AGROW>
	end
	
	fprintf('Func %d, best solver: %s, fval: %0.2E\n', ...
		i, bestSolvers, minErrs(i));
end
toc(comparisonStart);
end