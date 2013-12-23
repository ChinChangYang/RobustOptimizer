function y = genotype2phenotype(x)
% GENOTYPE2PHENOTYPE Transform genotype to phenotype for jadeeig
% See also PHENOTYPE2GENOTYPE, JADEEIGERR
[phenotype_lb, phenotype_ub] = jadeeig.chromosome.phenotypebounds;
y = phenotype_lb + (phenotype_ub - phenotype_lb) .* x;
end

