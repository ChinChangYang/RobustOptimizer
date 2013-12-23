function [X, Y, Z, C, h] = drawerrorcontour
%DRAWERRORCONTOUR Draw error contour
global solver measureOptions fitfun optimizedGenotype
measureOptions.D = 2;
measureOptions.Runs = 1;
measureOptions.MaxFunEvals = round(1e2);
measureOptions.LowerBounds = -100;
measureOptions.UpperBounds = 100;
solver = 'projadeeig';
fitfun = 'cec13_f12';
startTime = tic;
addprojectpath;
load optimprojadeeigD2M100;
optimizedGenotype = xmin;
[X, Y, Z, C, h] = drawcontour(@solverMeasure, 0, 0.01, 1);
saveFilename = sprintf('solver/soea/%s/%scontour.mat', solver, solver);
save(saveFilename, 'X', 'Y', 'Z');
toc(startTime);
end
function ret = solverMeasure(x)
global solver measureOptions fitfun optimizedGenotype
genotype = optimizedGenotype;
genotype(1) = x(1);
genotype(4) = x(2);
ret = computeerrors(fitfun, solver, genotype, measureOptions);
end