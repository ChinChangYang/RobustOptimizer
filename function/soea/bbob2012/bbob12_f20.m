function [Fval, Ftrue] = bbob12_f20(x)
% Schwefel with tridiagonal variable transformation
% last change: 08/12/17

persistent Fopt Xopt scales
persistent lastSize arrXopt arrScales arrSigns rseed

funcID = 20;
condition = 10;  % for linear transformation

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
% DIM-dependent initialization
if isempty(lastSize) || lastSize.DIM ~= DIM
	Xopt = 0.5 * sign(unif(DIM,rseed)'-0.5) * 4.2096874633;
	scales = sqrt(condition).^linspace(0, 1, DIM)';
end
% DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
	lastSize.POPSI = POPSI;
	lastSize.DIM = DIM;
	arrXopt = repmat(2*abs(Xopt), 1, POPSI);
	% arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI);
	arrScales = repmat(scales, 1, POPSI);
	arrSigns = repmat(sign(Xopt), 1, POPSI);  %
end

%----- TRANSFORMATION IN SEARCH SPACE -----
x = 2 * arrSigns .* x;  % makes the below boundary handling effective for coordinates
x(2:end,:) = x(2:end,:) + 0.25 * (x(1:end-1,:) - arrXopt(1:end-1,:));
x = 100 * (arrScales .* (x - arrXopt) + arrXopt);

%----- BOUNDARY HANDLING -----
xoutside = max(0, abs(x) - 500) .* sign(x);  % in [-500,500]
Fpen = 0.01 * sum(xoutside.^2, 1);  % penalty
Fadd = Fadd + Fpen;

%----- COMPUTATION core -----
% 0.01: values below one should be close to the global optimum
Ftrue = 0.01 * ((418.9828872724339) - mean(x.*sin(sqrt(abs(x))), 1));
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
	else  % if strcmpi(strinput, 'info')
		Ftrue = []; % benchmarkinfos(funcID);
		Fval = Fopt;
	end
end
end
