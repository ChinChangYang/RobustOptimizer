function optimize_de(D, fnum, maxfes)
fprintf('Optimizing CEC15_D%df%d...\n', D, fnum);
solver	= 'SPS_L_SHADE_EIG';
fitfun	= @(x) generic_solver_err(...
    x, ...
    solver, ...
    D, ...
    sprintf('cec15_f%d', fnum), ...
    D * 50);
ftarget = -1;
startTime = tic;

lb = [...
	0; ...		% NP - NPmin
	0.0; ...	% F
	0.0; ...	% CR
	0.0; ...	% ER
	0.0; ...	% p
	2; ...		% H
	1; ...		% Q
	1; ...		% Ar
	0.0; ...	% cw
	0.0; ...	% erw
	0.0; ...	% CRmin
	0.0; ...	% CRmax - CRmin
	4; ...		% NPmin
	0.0; ...	% crw
	eps; ...	% fw
	];

ub = [...
	20 * D; ...	% NP - NPmin
	1.0; ...	% F
	1.0; ...	% CR
	1.0; ...	% ER
	1.0; ...	% p
	20 * D; ...	% H
	20 * D; ...	% Q
	5; ...		% Ar
	1.0; ...	% cw
	1.0; ...	% erw
	1.0; ...	% CRmin
	1.0; ...	% CRmax - CRmin
	20 * D; ...	% NPmin
	1.0; ...	% crw
	1.0; ...	% fw
	];

prior1 = [...
	D * 19 - 4; ...	% NP - NPmin
	0.5; ...	% F
	0.5; ...	% CR
	1.0; ...	% ER
	0.11; ...	% p
	6; ...		% H
	64; ...		% Q
	2.6; ...	% Ar
	0.3; ...	% cw
	0.2; ...	% erw
	0.05; ...	% CRmin
	0.25; ...	% CRmax - CRmin
	4; ...		% NPmin
	0.1; ...	% crw
	0.1; ...	% fw
	];

prior2 = [...
	D * 18 - 4; ...	% NP - NPmin
	0.5; ...	% F
	0.5; ...	% CR
	0.0; ...	% ER
	0.11; ...	% p
	6; ...		% H
	20 * D; ... % Q
	2.6; ...	% Ar
	0.0; ...	% cw
	0.0; ...	% erw
	0.0; ...	% CRmin
	1.0; ...	% CRmax - CRmin
	4; ...		% NPmin
	0.1; ...	% crw
	0.1; ...	% fw
	];

NP = 2 * numel(prior1);
options.NP = NP;
options.initial.X = repmat(lb, 1, NP) + repmat(ub - lb, 1, NP) .* rand(numel(prior1), NP);
options.initial.X(:, 1) = prior1;
options.initial.X(:, 2) = prior2;
options.Q = 1;
options.ftarget = ftarget;
options.Display = 'iter';
options.Plotting = 'off';
options.EarlyStop = 'auto';
options.Noise	= true;

[xmin, fmin, out] = ...
	feval(solver, fitfun, lb, ub, maxfes, options); %#ok<ASGLU>

filename = sprintf('optimize_de_result_D%df%d.mat', D, fnum);
save(filename, 'xmin', 'fmin', 'out', 'solver', 'fitfun', 'lb', 'ub', ...
	'maxfes', 'options');

toc(startTime);
fprintf('Optimizing CEC15_D%df%d...done\n', D, fnum);
end