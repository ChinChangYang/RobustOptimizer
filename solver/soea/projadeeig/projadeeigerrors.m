function err = projadeeigerrors(x)
% PROJADEEIGERRORS Pro JADE/eig errors
% See also GENOTYPE2PHENOTYPE, PHENOTYPE2GENOTYPE
solver = 'projadeeig';

% Experiment settings
solverOptions = genotype2options(x, solver);
err = cheaperrors(solver, solverOptions);
end

