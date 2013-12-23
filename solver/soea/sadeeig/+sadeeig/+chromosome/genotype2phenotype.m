function y = genotype2phenotype(x)
% GENOTYPE2PHENOTYPE Transform genotype to phenotype for sadeeig
% See also PHENOTYPE2GENOTYPE, SADEEIGERROR
[phenotype_lb, phenotype_ub] = sadeeig.chromosome.phenotypebounds;
y = phenotype_lb + (phenotype_ub - phenotype_lb) .* x;
end

