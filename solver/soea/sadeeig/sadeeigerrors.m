function err = sadeeigerrors(x)
% SADEEIGERRORS SaDE/eig errors
% phenotype: [dimensionFactor, maxfunevalsFactor, R, LP]
% See also GENOTYPE2PHENOTYPE, PHENOTYPE2GENOTYPE
solver = 'sadeeig';

% Experiment settings
solverOptions = genotype2options(x, solver);
err = cheaperrors(solver, solverOptions);
end

