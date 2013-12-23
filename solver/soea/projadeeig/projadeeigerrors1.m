function ret = projadeeigerrors1(genotype)
%PROJADEEIGERRORS1 Error#1 of the Pro JADE/eig in cec13_f15
solver = 'projadeeig';
errorMeasure = 'cec13_f15';
ret = computeerrors(errorMeasure, solver, genotype);
end
