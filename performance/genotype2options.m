function solverOptions = genotype2options(x, solver)
% GENOTYPE2OPTOINS Transform genotype to options for a solver
% See also OPTIMIZESOLVER

[phenotype_lb, phenotype_ub] = ...
	eval(sprintf('%s.chromosome.phenotypebounds', solver));

D = numel(x);
x = reshape(x, 1, D);
y = phenotype_lb + (phenotype_ub - phenotype_lb) .* x; %#ok<NASGU>
solverOptions = ...
	eval(sprintf('%s.chromosome.phenotype2options(y)', solver));
end

