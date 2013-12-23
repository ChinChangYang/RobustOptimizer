function err = derand1eigerrors(x)
% DERAND1EIGERRORS DE/rand/1/eig errors
% See also GENOTYPE2PHENOTYPE, PHENOTYPE2GENOTYPE
solver = 'derand1eig';

% Experiment settings
solverOptions = genotype2options(x, solver);
err = cheaperrors(solver, solverOptions);
end

