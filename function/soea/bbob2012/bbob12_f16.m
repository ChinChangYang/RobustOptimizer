function [Fval, Ftrue] = bbob12_f16(x)
% Weierstrass, condition 100
% last change: 09/01/26
persistent Fopt Xopt scales linearTF rotation
persistent lastSize arrXopt rseed K aK bK F0

funcID = 16;
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
% DIM-dependent initialization
if isempty(lastSize) || lastSize.DIM ~= DIM
	Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation
	rotation = computeRotation(rseed+1e6, DIM);
	scales = (1/sqrt(condition)).^linspace(0, 1, DIM)';  % CAVE
	linearTF = diag(scales) * computeRotation(rseed, DIM);
	% decouple scaling from function definition
	linearTF = rotation * linearTF; % or computeRotation(rseed+1e3, DIM)
	K = 0:11; % number of summands, 20 in CEC2005, 10/12 saves 30% of time
	aK = 0.5.^K;
	bK = 3.^K;
	F0 = sum(aK .* cos(2*pi*bK * 0.5)); % optimal value
end
% DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
	lastSize.POPSI = POPSI;
	lastSize.DIM = DIM;
	arrXopt = repmat(Xopt, 1, POPSI);
	% arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI);
	% arrScales = repmat(scales, 1, POPSI);
end

%----- BOUNDARY HANDLING -----
xoutside = max(0, abs(x) - 5) .* sign(x);
Fpen = (10/DIM) * sum(xoutside.^2, 1);  % penalty
Fadd = Fadd + Fpen;

%----- TRANSFORMATION IN SEARCH SPACE -----
x = x - arrXopt;  % shift optimum to zero
x = rotation * x;  % no scaling here, because it would go to the arrExpo
x = monotoneTFosc(x);
x = linearTF * x;    % rotate/scale

%----- COMPUTATION core -----
Ftrue = zeros(1, POPSI);
for k = 1:POPSI
	% component-wise we have for x(i,k):
	%   sum(a.^K .* cos(2*pi* b.^K * (x(i,k)+0.5))) ...
	%     - DIM * sum(a.^K .* cos(2*pi*b.^K*0.5));
	Ftrue(k) = sum(cos(2*pi * (x(:,k)+0.5) * bK) * aK');  % the outer product *bK takes most time
end
Ftrue = 10 * (Ftrue/DIM - F0).^3;
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
