function options = phenotype2options(y)
%PHENOTYPE2OPTIONS Transform phenotype to options for jadeeig.

options.dimensionFactor = y(1);
options.R = y(2);
options.delta_CR = y(3);
options.delta_F = y(4);
options.p = y(5);
options.w = y(6);
options.TolFun = y(7);
options.FactorNP = y(8);
end
