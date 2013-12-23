function options = phenotype2options(y)
%PHENOTYPE2OPTIONS Transform phenotype to options for Pro JADE/eig.

options.dimensionFactor = y(1);
options.R = y(2);
options.delta_CR = y(3);
options.delta_F = y(4);
options.p = y(5);
options.w = y(6);
options.RdExp = y(7);
options.TolFun = y(8);
options.FactorNP = y(9);
end
