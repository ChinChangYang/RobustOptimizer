function x = phenotype2genotype( y )
% PHENOTYPE2GENOTYPE Transform phenotype to genotype for Pro JADE/eig
% See also GENOTYPE2PHENOTYPE, SADEEIGERROR
solver = 'projadeeig';
[phenotype_lb, phenotype_ub] = ...
	eval(sprintf('%s.chromosome.phenotypebounds', solver));
x = (y - phenotype_lb) ./ (phenotype_ub - phenotype_lb);
end
