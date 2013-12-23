function [Fval, Ftrue] = bbob12_f5(x)
% linear slope
% last change: 08/12/05
persistent Fopt Xopt scales
persistent lastSize arrXopt rseed

funcID = 5;
alpha = 100;

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
	Fopt =1* min(1000, max(-1000, (round(100*100*gauss(1,rseed)/gauss(1,rseed+1))/100)));
end
Fadd = Fopt;  % value to be added on the "raw" function value
% DIM-dependent initialization
if isempty(lastSize) || lastSize.DIM ~= DIM
	Xopt = 5 * sign(computeXopt(rseed, DIM));
	scales = - sign(Xopt) .* sqrt(alpha).^linspace(0, 1, DIM)';
end
Fadd = Fadd + 5 * sum(abs(scales));
% DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
	lastSize.POPSI = POPSI;
	lastSize.DIM = DIM;
	arrXopt = repmat(Xopt, 1, POPSI);
	% arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI);
	% arrScales = repmat(scales, 1, POPSI);
end

%----- BOUNDARY HANDLING -----
% move "too" good coordinates back into domain
idx_out_of_bounds = x .* arrXopt > 25;  % 25 == 5 * 5
x(idx_out_of_bounds) = sign(x(idx_out_of_bounds)) * 5;

%----- TRANSFORMATION IN SEARCH SPACE -----

%----- COMPUTATION core -----
Ftrue = scales' * x;
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
	else  % if strcmpi(strinput, 'info')
		Ftrue = []; % benchmarkinfos(funcID);
		Fval = Fopt;
	end
end
end
