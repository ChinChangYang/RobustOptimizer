function [Fval, Ftrue] = bbob12_f7(x)
% step-ellipsoid, condition 100
% last change: 08/12/05
persistent Fopt Xopt scales linearTF rotation
persistent lastSize arrXopt rseed

funcID = 7;
condition = 100;
alpha = 10;  % inner rounding procedure
rrseed = funcID;

%----- CHECK INPUT -----
if ischar(x) % return Fopt Xopt or linearTF on string argument
	flginputischar = 1;
	strinput = x;
	if nargin < 2
		DIM = 2;
	end
	x = ones(DIM,1);  % setting all persistent variables
else
	flginputischar = 0;
end
% from here on x is assumed a numeric variable
[DIM, POPSI] = size(x);  % dimension, pop-size (number of solution vectors)
if DIM == 1
	error('1-D input is not supported');
end

%----- INITIALIZATION -----
if nargin > 2     % set seed depending on trial index
	Fopt = [];      % clear previous settings for Fopt
	lastSize = [];  % clear other previous settings
	rseed = rrseed + 1e4 * ntrial;
elseif isempty(rseed)
	rseed = rrseed;
end
if isempty(Fopt)
	Fopt =1* min(1000, max(-1000, (round(100*100*gauss(1,rseed)/gauss(1,rseed+1))/100)));
end
Fadd = Fopt;  % value to be added on the "raw" function value
% DIM-dependent initialization
if isempty(lastSize) || lastSize.DIM ~= DIM
	Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation
	rotation = computeRotation(rseed+1e6, DIM);
	scales = condition.^linspace(0, 1, DIM)';
	linearTF = diag(sqrt(condition/10).^linspace(0, 1, DIM)) * computeRotation(rseed, DIM);
end
% DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
	lastSize.POPSI = POPSI;
	lastSize.DIM = DIM;
	arrXopt = repmat(Xopt, 1, POPSI);
end

%----- BOUNDARY HANDLING -----
xoutside = max(0, abs(x) - 5) .* sign(x);
Fpen = 1 * sum(xoutside.^2, 1);  % penalty
Fadd = Fadd + Fpen;

%----- TRANSFORMATION IN SEARCH SPACE -----
x = x - arrXopt;  % shift optimum to zero
x = linearTF * x;  % rotate while assuming that Xopt == 0
x1 = x(1,:);
idx = abs(x) > 0.5;
x(idx) = round(x(idx));
x(~idx) = round(alpha*x(~idx))/alpha;
x = rotation * x;  % rotate while assuming that Xopt == 0

%----- COMPUTATION core -----
Ftrue = 0.1 * max(1e-4*abs(x1), scales' * x.^2);
% Ftrue(Ftrue == 0) = 1e-4 * x1(Ftrue == 0).^2;
Fval = Ftrue;  % without noise

%----- FINALIZE -----
Ftrue = Ftrue + Fadd;
Fval = Fval + Fadd;

%----- RETURN INFO -----
if flginputischar
	if strcmpi(strinput, 'xopt')
		Fval = Fopt;
		Ftrue = Xopt;
	elseif strcmpi(strinput, 'linearTF')
		Fval = Fopt;
		Ftrue = {};
		Ftrue{1} = linearTF;
		Ftrue{2} = rotation;
	else  % if strcmpi(strinput, 'info')
		Ftrue = []; % benchmarkinfos(funcID);
		Fval = Fopt;
	end
end
end
