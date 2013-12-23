function illustrate_clop_functions
%ILLUSTRATE_CLOP_FUNCTIONS Illustrate CLOP functions
fitfuns = {'cloplog', 'clopflat', 'cloppower', 'clopangle', 'clopstep'};

for i = 1 : numel(fitfuns)
	options.fitfun = fitfuns{i};
	options.T = 200;	
	fitfun = @(x)functionsampler(x, options);
	figure;
	plotfunction(fitfun, -1, 0.01, 1);
end

fitfuns = {'cloplognoisefree', 'cloplog', ...
	'clopflatnoisefree', 'clopflat', ...
	'cloppowernoisefree', 'cloppower', ...
	'clopanglenoisefree', 'clopangle', ...
	'clopstepnoisefree', 'clopstep', ...
	'cloprosenbrocknoisefree', 'cloprosenbrock', ...
	'clopcorrelatednoisefree', 'clopcorrelated'};

for i = 1 : numel(fitfuns)
	options.fitfun = fitfuns{i};
	options.T = 2;
	fitfun = @(x)functionsampler(x, options);
	figure;
	drawcontour(fitfun, -1, 0.02, 1);
end
end

