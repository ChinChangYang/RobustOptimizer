function [Fval, Ftrue] = bbob12_f21(x)
% Gallagher with 101 Gaussian peaks, condition up to 1000, one global rotation
% CAVE: the difficulty might significantly depend on the
%       random realization for the location of optima.
% last change: 09/01/26
persistent Fopt Xopt rotation Xlocal peakvalues
persistent lastSize arrScales rseed

funcID = 21;
maxcondition = 1000;  % for linear transformation
fitvalues = [1.1, 9.1];    % range for raw values, 10 is optimal
nhighpeaks = 101;

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
	rotation = computeRotation(rseed, DIM);
	% compute scaling for each optimum
	arrCondition = maxcondition.^linspace(0,1,nhighpeaks-1);
	[~, idx] = sort(unif(nhighpeaks-1, rseed));  % random permutation
	arrCondition = [sqrt(maxcondition) arrCondition(idx)];
	arrScales = zeros(nhighpeaks, DIM);
	for i = 1:nhighpeaks  % generation of cov matrices
		s = arrCondition(i).^linspace(-0.5, 0.5, DIM);
		[~, idx] = sort(unif(DIM, rseed+1e3*(i-1))); % permutation instead of rotation
		arrScales(i,:) = s(idx); % this is inverse Cov
	end
	% compute peak values, 10 is global optimum
	peakvalues = [10 linspace(fitvalues(1), fitvalues(2), nhighpeaks-1)];
end
% DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
	lastSize.POPSI = POPSI;
	lastSize.DIM = DIM;
	Xlocal = rotation * reshape(10*unif(DIM*nhighpeaks, rseed)-5, ...
		DIM, nhighpeaks);
	% global optimum not too close to boundary
	Xlocal(:,1) =1* 0.8 * Xlocal(:,1);
	Xopt = rotation' * Xlocal(:,1);
end

%----- BOUNDARY HANDLING -----
xoutside = max(0, abs(x) - 5) .* sign(x);
Fpen = 1 * sum(xoutside.^2, 1);  % penalty
Fadd = Fadd + Fpen;

%----- TRANSFORMATION IN SEARCH SPACE -----
x = rotation * x;

%----- COMPUTATION core -----
fac = -0.5/DIM; % ??xopt + alpha*ones versus alpha becomes invariant like this
f = NaN(nhighpeaks,POPSI);
if POPSI < 0.5*nhighpeaks
	for k = 1:POPSI
		xx = repmat(x(:,k), 1, nhighpeaks) - Xlocal;
		f(:,k) = peakvalues' .* exp(fac*(sum(arrScales'.*xx.^2, 1)))';
	end
else
	for i = 1:nhighpeaks  % CAVE: POPSI = 1e4 gets out of memory
		xx = (x - repmat(Xlocal(:,i), 1, POPSI)); % repmat could be done once
		f(i,:) = peakvalues(i) * exp(fac*(arrScales(i,:)*xx.^2));
	end
end
% f is in [0,10], 10 is best
Ftrue = monotoneTFosc(10 - max(f, [], 1)).^2;
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
		Ftrue{2} = rotation;
	else  % if strcmpi(strinput, 'info')
		Ftrue = []; % benchmarkinfos(funcID);
		Fval = Fopt;
	end
end
end
