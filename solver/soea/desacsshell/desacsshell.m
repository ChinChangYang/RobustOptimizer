function [xmin, fmin, out] = desacsshell(fitfun, lb, ub, maxfunevals, options)
% DESACSSHELL Shell of DESaCS
% DESACSSHELL(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DESACSSHELL(..., options) minimize the function by solver options.

x0 = (ub + lb) / 2;
insteps = ub - lb;
defaultOptions = desacsdefoptions;
options.maxfunevals = maxfunevals;
options = setdefoptions(options, defaultOptions);

if isfield(options, 'Display')
	if strcmp(options.Display, 'iter')
		options.plotting = true;
		options.de_more_plot = true;
		options.contourl = -100;
		options.contourr = 5;
		options.contouru = 100;
	end
end

[~, ~, ~, ~, out, bestever] ...
	= desacs(fitfun, x0, insteps, options);

xmin = bestever.xmin;
fmin = bestever.fmin;
out.bestever.xmin = xmin;
out.bestever.fmin = fmin;
end
