function [Fval, Ftrue] = bbob12_f24(x)
% Lunacek bi-Rastrigin, condition 100
% in PPSN 2008, Rastrigin part rotated and scaled
% last change: 08/12/17
persistent Fopt Xopt scales linearTF rotation
persistent lastSize arrScales rseed

funcID = 24;
condition = 100;  % for linear transformation

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
	rseed = funcID + 1e4 * ntrial;
elseif isempty(rseed)
	rseed = funcID;
end
if isempty(Fopt)
	Fopt =1* min(1000, max(-1000, round(100*100*gauss(1,rseed)/gauss(1,rseed+1))/100));
end
Fadd = Fopt;  % value to be added on the "raw" function value
mu1 = 2.5;  % optimum shift
% DIM-dependent initialization
if isempty(lastSize) || lastSize.DIM ~= DIM
	Xopt = 0.5 * mu1 * sign(gauss(DIM, rseed))' .* ones(DIM, 1);
	rotation = computeRotation(rseed+1e6, DIM);
	scales = sqrt(condition).^linspace(0, 1, DIM)';
	linearTF = diag(scales) * computeRotation(rseed, DIM);
	% decouple scaling from function definition
	linearTF = rotation * linearTF; % or computeRotation(rseed+1e3, DIM)
end
% DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
	lastSize.POPSI = POPSI;
	lastSize.DIM = DIM;
	arrScales = repmat(2 * sign(Xopt), 1, POPSI);  % makes up for Xopt
end

%----- BOUNDARY HANDLING -----
xoutside = max(0, abs(x) - 5) .* sign(x);
Fpen = 1e4 * sum(xoutside.^2, 1);  % penalty
Fadd = Fadd + Fpen;

%----- TRANSFORMATION IN SEARCH SPACE -----
x = arrScales .* x;

%----- COMPUTATION core -----
s = 1 - 0.5/(sqrt(DIM+20)-4.1);  % tested up to DIM = 160 p in [0.25,0.33]
d = 1;  % shift [1,3], smaller is more difficult
mu2 = -sqrt((mu1^2 - d) / s);  % TODO: re-check this
Ftrue = min(sum((x-mu1).^2, 1), d*DIM + s * sum((x-mu2).^2, 1));
Ftrue = Ftrue + 10*(DIM - sum(cos(2*pi*(linearTF * (x-mu1))), 1));
Fval = Ftrue;  % without noise

%----- NOISE -----

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
