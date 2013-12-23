function options = phenotype2options(y)
% PHENOTYPE2OPTIONS Transform phenotype to options for derand1eig.
% phenotype: [dimensionFactor, maxfunevalsFactor, R, CR, F]
options.dimensionFactor = y(1);
options.maxfunevalsFactor = y(2);
options.R = y(3);
options.CR = y(4);
options.F = y(5);
end
