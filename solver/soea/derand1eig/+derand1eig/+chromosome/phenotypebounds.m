function [phenotype_lb, phenotype_ub] = phenotypebounds
%PHENOTYPEBOUNDS Get lower and upper bounds of the phenotype for derand1eig
phenotype_lb = [1, 0.0005, 0.0, 0.05, 0.4];
phenotype_ub = [10, 0.05, 1.0, 0.95, 1.0];
end

