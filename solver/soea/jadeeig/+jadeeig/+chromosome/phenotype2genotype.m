function x = phenotype2genotype( y )
% PHENOTYPE2GENOTYPE Transform phenotype to genotype for jadeeig
% See also GENOTYPE2PHENOTYPE, JADEEIGERR
[phenotype_lb, phenotype_ub] = jadeeig.chromosome.phenotypebounds;
x = (y - phenotype_lb) ./ (phenotype_ub - phenotype_lb);
end

