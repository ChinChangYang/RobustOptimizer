function y = genotype2phenotype(x)
% GENOTYPE2PHENOTYPE Transform genotype to phenotype for Pro JADE/eig
% See also PHENOTYPE2GENOTYPE, SADEEIGERROR
solver = 'projadeeig';
[phenotype_lb, phenotype_ub] = ...
	eval(sprintf('%s.chromosome.phenotypebounds', solver));
y = phenotype_lb + (phenotype_ub - phenotype_lb) .* x;
end

