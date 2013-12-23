function [Fval, Ftrue] = bbob12_f18(x)
% Schaffers F7 with asymmetric non-linear transformation, condition 1000
% last change: 08/12/17
persistent Fopt Xopt scales linearTF rotation
persistent lastSize arrXopt arrExpo rseed

condition = 1e3;  % for linear transformation
beta = 0.5;      %
rrseed = 17;

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
	Fopt =1* min(1000, max(-1000, round(100*100*gauss(1,rseed)/gauss(1,rseed+1))/100));
end
Fadd = Fopt;  % value to be added on the "raw" function value
% DIM-dependent initialization
if isempty(lastSize) || lastSize.DIM ~= DIM
	Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation
	rotation = computeRotation(rseed+1e6, DIM);
	scales = sqrt(condition).^linspace(0, 1, DIM)';
	linearTF = diag(scales) * computeRotation(rseed, DIM);
	% decouple scaling from function definition
	% linearTF = rotation * linearTF; % or computeRotation(rseed+1e3, DIM)
end
% DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
	lastSize.POPSI = POPSI;
	lastSize.DIM = DIM;
	arrXopt = repmat(Xopt, 1, POPSI);
	arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI);
	% arrScales = repmat(scales, 1, POPSI);
end

%----- BOUNDARY HANDLING -----
xoutside = max(0, abs(x) - 5) .* sign(x);
Fpen = 1e1 * sum(xoutside.^2, 1);  % penalty
Fadd = Fadd + Fpen;

%----- TRANSFORMATION IN SEARCH SPACE -----
x = x - arrXopt;  % shift optimum to zero
x = rotation * x;  % no scaling here, because it would go to the arrExpo
% x = monotoneTFosc(x);
idx = x > 0;
x(idx) = x(idx).^(1 + arrExpo(idx) .* sqrt(x(idx)));  % smooth in zero
% x = rotation' * x;   % back-rotation
% x = arrScales .* x;  % scaling, Xopt should be 0
x = linearTF * x;    % rotate/scale

%----- COMPUTATION core -----
s = x(1:DIM-1,:).^2 + x(2:DIM,:).^2;
Ftrue = mean(s.^0.25 .* (sin(50*s.^0.1).^2+1), 1).^2; % ^2 for ftarget
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
