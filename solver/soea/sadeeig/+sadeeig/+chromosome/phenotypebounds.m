function [phenotype_lb, phenotype_ub] = phenotypebounds
%PHENOTYPEBOUNDS Get lower and upper bounds of the phenotype for sadeeig
phenotype_lb = [1, 1e-5, 0, 5];
phenotype_ub = [10, 1e-3, 1.0, 100];
end
