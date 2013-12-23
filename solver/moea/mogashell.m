function [xmin, fmin] = mogashell(fitfun, lb, ub, maxfunevals, options)
% MOGASHELL Shell of multiobjective GA
fhandle = str2func(fitfun);
D = numel(lb);
NP = max(2 * D, floor(maxfunevals / 3e2));
Gmax = floor(maxfunevals / NP);

LHS = lhsdesign(NP, D, 'iteration', 50);
Pop = zeros(NP, D);
for i = 1 : NP
	Pop(i, :) = lb' + (ub - lb)' .* LHS(i, :);
end

if nargin == 5
	if isstruct(options) && isfield(options, 'InitialPopulation')
		Pop(end, :) = options.InitialPopulation;
	end
end

options = gaoptimset(@gamultiobj);
options = gaoptimset(options, 'Display', 'iter', ...
	'PopulationSize', NP, ...
	'Generations', Gmax, ...
	'InitialPopulation', Pop, ...
	'PopInitRange', [lb'; ub'], ...
	'TolFun', 0, ...
	'StallGenLimit', Inf);

[xmin, fmin] = gamultiobj(fhandle, D, [], [], [], [], lb, ub, options);

end

