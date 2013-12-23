function options = getoptimizedoptions(solver, iOptions, D, M)
%GETOPTIMIZEDOPTIONS Get optimized options

addprojectpath;

if nargin <= 2
	datafilename = sprintf('optim%s.mat', solver);
elseif nargin == 4
	datafilename = sprintf('optim%sD%dM%d.mat', solver, D, M);
end

optimData = load(datafilename);
xmin = optimData.xmin; %#ok<NASGU>

exp = sprintf('%s.chromosome.genotype2phenotype(xmin(%d, :))', ...
	solver, iOptions);
phenotype = eval(exp); %#ok<NASGU>
exp = sprintf('%s.chromosome.phenotype2options(phenotype)', solver);
options = eval(exp);
end

