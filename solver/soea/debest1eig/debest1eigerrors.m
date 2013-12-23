function err = debest1eigerrors(x)
% DEBEST1EIGERRORS DE/best/1/eig errors
% See also GENOTYPE2PHENOTYPE, PHENOTYPE2GENOTYPE
solver = 'debest1eig';

% Experiment settings
solverOptions = genotype2options(x, solver);
err = cheaperrors(solver, solverOptions);
end

