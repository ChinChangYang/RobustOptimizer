function options = phenotype2options(y)
%PHENOTYPE2OPTIONS Transform phenotype to options for sadeeig.
options.dimensionFactor = y(1);
options.maxfunevalsFactor = y(2);
options.R = y(3);
options.LP = y(4);
end
