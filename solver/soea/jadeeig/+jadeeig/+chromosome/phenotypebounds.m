function [phenotype_lb, phenotype_ub] = phenotypebounds
%PHENOTYPEBOUNDS Get lower and upper bounds of the phenotype for jadeeig
phenotype_lb = [1, 0, 0.01, 0.01, 0.05, 0.05, 1e-10, 1];
phenotype_ub = [4, 1, 0.3, 0.3, 0.5, 0.5, 1e-8, 3];
end

