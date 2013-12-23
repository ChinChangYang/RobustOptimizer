function solutionerrordistribution
%SOLUTIONERRORDISTRIBUTION Show distribution of solution errors
global solver measureOptions fitfun
solver = 'projadeeig';
measureOptions.D = 2;
measureOptions.Runs = 1;
measureOptions.MaxFunEvals = round(1e2);
measureOptions.LowerBounds = -100;
measureOptions.UpperBounds = 100;
fitfun = 'cec13_f12';

startTime = tic;
addprojectpath;
load optimprojadeeigD2M100;
T = round(1e6 / measureOptions.MaxFunEvals);
solutionError = zeros(1, T);
for t = 1 : T
	solutionError(t) = solverMeasure(xmin);
	hist(solutionError(1:t), 100);
	fprintf('Mean of solution errors: %0.4E, Std. of solution errors: %0.4E\n', ...
		mean(solutionError(1:t)), std(solutionError(1:t)));
	pause(0.01);
end
xlabel('Solution error');
ylabel('Frequency');
sTitle = ...
	sprintf(['Histogram of solution errors of Pro JADE/eig in function' ...
	' 12 for D=%d, M=%d, and T=%d'], ...
	measureOptions.D, measureOptions.MaxFunEvals, T);
	
title(sTitle);
toc(startTime);
end
function ret = solverMeasure(genotype)
global solver measureOptions fitfun
ret = computeerrors(fitfun, solver, genotype, measureOptions);
end
