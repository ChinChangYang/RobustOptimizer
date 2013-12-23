function ret = projadeeigerrors2(genotype)
%PROJADEEIGERRORS2 Error#2 of the Pro JADE/eig in cec13_f23
solver = 'projadeeig';
errorMeasure = 'cec13_f23';
ret = computeerrors(errorMeasure, solver, genotype);
end

