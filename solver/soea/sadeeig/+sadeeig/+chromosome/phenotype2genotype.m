function x = phenotype2genotype( y )
% PHENOTYPE2GENOTYPE Transform phenotype to genotype for sadeeig
% See also GENOTYPE2PHENOTYPE, SADEEIGERROR
[phenotype_lb, phenotype_ub] = sadeeig.chromosome.phenotypebounds;
x = (y - phenotype_lb) ./ (phenotype_ub - phenotype_lb);
end
