function x = phenotype2genotype( y )
% PHENOTYPE2GENOTYPE Transform phenotype to genotype for derand1eig
% See also GENOTYPE2PHENOTYPE, DERAND1EIGERRORS
solver = 'derand1eig';
[phenotype_lb, phenotype_ub] = ...
	eval(sprintf('%s.chromosome.phenotypebounds', solver));
x = (y - phenotype_lb) ./ (phenotype_ub - phenotype_lb);
end
