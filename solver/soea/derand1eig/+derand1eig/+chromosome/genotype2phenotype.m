function y = genotype2phenotype(x)
% GENOTYPE2PHENOTYPE Transform genotype to phenotype for derand1eig
% See also PHENOTYPE2GENOTYPE, DERAND1EIGERRORS
solver = 'derand1eig';
[phenotype_lb, phenotype_ub] = ...
	eval(sprintf('%s.chromosome.phenotypebounds', solver));
y = phenotype_lb + (phenotype_ub - phenotype_lb) .* x;
end

