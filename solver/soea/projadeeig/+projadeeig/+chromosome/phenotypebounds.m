function [phenotype_lb, phenotype_ub] = phenotypebounds
%PHENOTYPEBOUNDS Get lower and upper bounds of the phenotype for Pro
%JADE/eig
phenotype_lb = [1, 0, 0.01, 0.01, 0.01, 0.01, 0, 1e-11, 1];
phenotype_ub = [10, 1, 0.5, 0.5, 0.5, 0.5, 3, 1e-8, 4];
end
