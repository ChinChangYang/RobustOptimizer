% VAL = BENCHMARKS(X, FUNCID)
% VAL = BENCHMARKS(X, STRFUNC)
%    Input: 
%       X -- solution column vector or matrix of column vectors
%       FUNCID -- number of function to be executed with X as input,
%                 by default 8. 
%       STRFUNC -- function as string to be executed with X as input
%    Output: function value(s) of solution(s)
%    Examples: 
%      F = BENCHMARKS([1 2 3]', 17); 
%      F = BENCHMARKS([1 2 3]', 'f1'); 
% 
% NBS = BENCHMARKS() 
% NBS = BENCHMARKS('FunctionIndices') 
%    Output: 
%      NBS -- array of valid benchmark function numbers, 
%             presumably 1:24
%
% FHS = BENCHMARKS('handles')
%    Output: 
%      FHS -- cell array of function handles
%    Examples:
%      FHS = BENCHMARKS('handles');  
%      f = FHS{1}(x);  % evaluates x on the sphere function f1
%      f = feval(FHS{1}, x);  % ditto 
%
% see also: functions FGENERIC, BENCHMARKINFOS, BENCHMARKSNOISY

% Authors (copyright 2009): Nikolaus Hansen, Raymond Ros, Steffen Finck
%    Version = 'Revision: $Revision: 3162 $'
%    Last Modified: $Date: 2009-02-09 19:22:42 +0100 (Mon, 09 Feb 2009) $

% INTERFACE OF BENCHMARK FUNCTIONS
% FHS = BENCHMARKS('handles');
% FUNC = FHS{1};
%
% [FVALUE, FTRUE] = FUNC(X)
% [FVALUE, FTRUE] = FUNC(X, [], IINSTANCE)
%   Input: X -- matrix of column vectors
%          IINSTANCE -- instance number of the function, sets function
%             instance (XOPT, FOPT, rotation matrices,...)
%             up until a new number is set, or the function is
%             cleared. Default is zero.
%   Output: row vectors with function value for each input column
%     FVALUE -- function value
%     FTRUE -- noise-less, deterministic function value
% [FOPT STRFUNCTION] = FUNC('any_even_empty_string', ...)
%   Output:
%     FOPT -- function value at optimum
%     STRFUNCTION -- not yet implemented: function description string, ID before first whitespace
% [FOPT STRFUNCTION] = FUNC('any_even_empty_string', DIM, NTRIAL)
%   Sets rotation matrices and xopt depending on NTRIAL (by changing the random seed). 
%   Output:
%     FOPT -- function value at optimum
%     STRFUNCTION -- not yet implemented: function description string, ID before first whitespace
% [FOPT, XOPT] = FUNC('xopt', DIM)
%   Output:
%     FOPT -- function value at optimum XOPT
%     XOPT -- optimal solution vector in DIM-D
% [FOPT, MATRIX] = FUNC('linearTF', DIM)  % might vanish in future
%   Output:
%     FOPT -- function value at optimum XOPT
%     MATRIX -- used transformation matrix 


%%%-------------------------------------------------------------%%%
function [res, res2] = benchmarks(x, strfun, DIM) 
%
  Nfcts = 24;

  % return valid function IDs (ie numbers)
  if nargin < 1 || ( ...
      ischar(x) && (strcmpi(x, 'noisefreeFunctionIndices') || strcmpi(x, 'FunctionIndices')))
    res = [];
    for i = 1:100
      try  % exist does not work under Octave
        eval(['f' num2str(i) '([1; 2]);']);
        res(end+1) = i;
      catch
        if i < Nfcts
          disp(sprintf('execution of function %d produces an error', i));
        end
      end
    end 
    if length(res) > Nfcts
      disp(res);
      error([num2str(length(res)) ' > 24 functions f1, f2,...,f24' ...
             ' are defined. Try to clear the workspace']);
    elseif length(res) < Nfcts
      disp(res);
      error([num2str(length(res)) ' < 24 functions are defined. ' ...
             'There might be an execution error ']);
    end
  elseif isnumeric(x)
    if nargin < 2
      res = f8(x); 
      return;
    end
    if isnumeric(strfun)
      res = feval(['f' num2str(strfun)], x);
    else
      res = feval(strfun, x);
    end
  % return function handles
  elseif ischar(x) && strcmpi(x, 'handles')
    res = {}; res2 = []; 
    for i = 1:100
      try  % exist does not work under Octave
        eval(['f' num2str(i) '([1; 2]);']);
        res{i} = str2func(['f' num2str(i)]); % eval(['@f' num2str(i)]);
      catch
        if i < Nfcts 
          disp(sprintf('execution of function %d produces an error', i));
        end
      end
    end
    % res2 = benchmarkinfos();
    return;
  elseif ischar(x) && strcmpi(x, 'test')
    error('to be implemented');
    % TODO
    for DIM = [2,3,5,10,20,40]
      % generate test points also outside the domain 
      % pop = generateTestPoints
      for ifun = benchmarks('FunctionIndices')
        % evaluate all functions at all the points, using fgeneric
      end
    end
    % dump the result 
  elseif ischar(x) && strcmpi(x, 'fopt')
    if isnumeric(strfun)
      res = feval(['f' num2str(strfun)], x);
    else
      res = feval(strfun, x);
    end
  elseif ischar(x) && strncmpi(x, 'i', 1)
    res = benchmarkinfos(strfun); 
    if nargin > 1
      res2 = res;
      res = res2{strfun};
    end
  else
    error('not a valid case of input');
  end
end

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f1(x, DIM, ntrial)
% sphere function 
% last change: 08/12/02
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrExpo rseed

  funcID = 1; 

  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    % rotation = computeRotation(rseed+1e6, DIM); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2,1); 
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f2(x, DIM, ntrial)
% separable ellipsoid with monotone transformation, condition 1e6
% last change: 08/12/02
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrExpo rseed

  funcID = 2; 
  condition = 1e6;  
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    scales = condition.^linspace(0, 1, DIM)'; 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = scales' * monotoneTFosc(x).^2; 
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

end % function


%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f3(x, DIM, ntrial)
% Rastrigin with monotone transformation separable "condition" 10
% last change: 08/12/05
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 3; 
  condition = 10; % for linear transformation
  beta = 0.2; 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    scales = sqrt(condition).^linspace(0, 1, DIM)'; 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrScales = repmat(scales, 1, POPSI); 
    arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI); 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = monotoneTFosc(x);
  idx = x > 0;
  x(idx) = x(idx).^(1 + arrExpo(idx) .* sqrt(x(idx)));  % smooth also in zero
  x = arrScales .* x;  
  
  %----- COMPUTATION core -----
  Ftrue = 10 * (DIM - sum(cos(2*pi*x), 1)) + sum(x.^2, 1); 
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f4(x, DIM, ntrial)
% skew Rastrigin-Bueche, condition 10, skew-"condition" 100
% last change: 08/12/05
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 4; 
  condition = 10; % for linear transformation
  alpha = 100;
  maxindex = inf; % 1:2:min(DIM,maxindex) are the skew variables
  rrseed = 3; 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt(1:2:min(DIM,maxindex)) = abs(Xopt(1:2:min(DIM,maxindex)));
    scales = sqrt(condition).^linspace(0, 1, DIM)'; 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
    % arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI); 
    arrScales = repmat(scales, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 1e2 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = monotoneTFosc(x); 
  idx = false(DIM, POPSI); 
  idx(1:2:min(DIM,maxindex), :) = x(1:2:min(DIM,maxindex), :) > 0; 
  x(idx) = sqrt(alpha)*x(idx);
  x = arrScales .* x;  % scale while assuming that Xopt == 0 

  %----- COMPUTATION core -----
  Ftrue = 10 * (DIM - sum(cos(2*pi*x), 1)) + sum(x.^2, 1);
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f5(x, DIM, ntrial)
% linear slope 
% last change: 08/12/05
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 5; 
  alpha = 100;  
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = []; % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f6(x, DIM, ntrial)
% attractive sector function
% last change: 08/12/17
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 6; 
  condition = 10;  % for linear transformation
  alpha = 100;      % 
  rrseed = funcID; 

  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    linearTF = rotation * linearTF; % or computeRotation(rseed+1e3, DIM)
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = linearTF * x;    % rotate/scale

  %----- COMPUTATION core -----
  idx = x.*arrXopt > 0;
  x(idx) = alpha*x(idx);
  Ftrue = monotoneTFosc(sum(x.^2,1)).^0.9; 
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f7(x, DIM, ntrial)
% step-ellipsoid, condition 100
% last change: 08/12/05
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 7; 
  condition = 100; 
  alpha = 10;  % inner rounding procedure
  rrseed = funcID; 

  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f8(x, DIM, ntrial)
% Rosenbrock, non-rotated
% last change: 08/12/10
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 8; 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* 0.75 * computeXopt(rseed, DIM); 
    scales = max(1, sqrt(DIM) / 8.);  % only needed for rotated case
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = scales * x;  % rotate/scale while assuming that Xopt == 0 
  x = x + 1;  % shift zero to factual optimum 1

  %----- COMPUTATION core -----
  Ftrue = 1e2 * sum((x(1:end-1, :).^2 - x(2:end, :)).^2, 1) ...
          + sum((x(1:end-1, :) - 1).^2, 1);
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

end % function


%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f9(x, DIM, ntrial)
% Rosenbrock, rotated
% last change: 08/12/10
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 9; 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    scale = max(1, sqrt(DIM) / 8.); 
    linearTF = scale * computeRotation(rseed, DIM); 
    Xopt = linearTF' * 0.5*ones(DIM,1) / scale^2;
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = linearTF * x + 0.5; 

  %----- COMPUTATION core -----
  Ftrue = 1e2*sum((x(1:end-1,:).^2 - x(2:end,:)).^2,1) + sum((x(1:end-1,:)-1).^2,1);
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f10(x, DIM, ntrial)
% ellipsoid with monotone transformation, condition 1e6
% last change: 08/12/10
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 10; 
  condition = 1e6;  % for linear transformation
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    rotation = computeRotation(rseed+1e6, DIM); 
    scales = condition.^linspace(0, 1, DIM)'; 
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

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;  % no scaling here, because it would go to the arrExpo
  x = monotoneTFosc(x);

  %----- COMPUTATION core -----
  Ftrue = scales' * x.^2; 
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f11(x, DIM, ntrial)
% discus (tablet) with monotone transformation, condition 1e6
% last change: 08/12/10
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 11; 
  condition = 1e6;

  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    rotation = computeRotation(rseed+1e6, DIM); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;  % no scaling here, because it would go to the arrExpo
  x = monotoneTFosc(x); 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2,1) + (condition-1) * x(1,:).^2;
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f12(x, DIM, ntrial)
% bent cigar with asymmetric space distortion, condition 1e6 
% last change: 08/12/11
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 12; 
  condition = 1e6; 
  beta = 0.5;  % distortion of space
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed+1e6, DIM); % function ID is seed for rotation 
    rotation = computeRotation(rseed+1e6, DIM); 
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

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;  % no scaling here, because it would go to the arrExpo
  idx = x > 0;
  x(idx) = x(idx).^(1 + arrExpo(idx) .* sqrt(x(idx)));  % smooth in zero
  x = rotation * x;  

  %----- COMPUTATION core -----
  Ftrue = condition * sum(x.^2, 1) + (1-condition) * x(1,:).^2;
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f13(x, DIM, ntrial)
% sharp ridge
% last change: 08/12/11
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 13; 
  condition = 10; 
  alpha = 100;  % slope 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    rotation = computeRotation(rseed+1e6, DIM); 
    scales = sqrt(condition).^linspace(0, 1, DIM)'; 
    linearTF = diag(scales) * computeRotation(rseed, DIM); 
    % decouple scaling from function definition
    linearTF = rotation * linearTF; 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = linearTF * x;    % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = x(1,:).^2 + alpha * sqrt(sum(x(2:end,:).^2, 1));
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f14(x, DIM, ntrial)%
% sum of different powers, between x^2 and x^6
% last change: 08/12/11
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 14; 
  alpha = 4; 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
    rotation = computeRotation(rseed+1e6, DIM); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
    arrExpo = repmat(2 + alpha * linspace(0, 1, DIM)', 1, POPSI); 
    % arrScales = repmat(scales, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;    % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = sqrt(sum(abs(x).^arrExpo, 1)); 
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f15(x, DIM, ntrial)
% Rastrigin with asymmetric non-linear distortion, "condition" 10
% last change: 08/12/26
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 15; 
  condition = 10;  % for linear transformation
  beta = 0.2;  
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    Xopt =1* computeXopt(rseed, DIM); % function ID is seed for rotation 
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
    arrXopt = repmat(Xopt, 1, POPSI); 
    arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;  % no scaling here, because it would go to the arrExpo
  x = monotoneTFosc(x); 
  idx = x > 0;
  x(idx) = x(idx).^(1 + arrExpo(idx) .* sqrt(x(idx)));  % smooth in zero
  x = linearTF * x;     % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = 10 * (DIM - sum(cos(2*pi*x), 1)) + sum(x.^2, 1); 
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f16(x, DIM, ntrial)
% Weierstrass, condition 100
% last change: 09/01/26
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed K aK bK F0

  funcID = 16; 
  condition = 100;  % for linear transformation
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    K = [0:11]; % number of summands, 20 in CEC2005, 10/12 saves 30% of time
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f17(x, DIM, ntrial)
% Schaffers F7 with asymmetric non-linear transformation, condition 10
% last change: 08/12/17
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 17; 
  condition = 10;  % for linear transformation
  beta = 0.5;      % 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
  Fpen = 1e1 * sum(xoutside.^2, 1);  % parallel to f18
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f18(x, DIM, ntrial)
% Schaffers F7 with asymmetric non-linear transformation, condition 1000
% last change: 08/12/17
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 18; 
  condition = 1e3;  % for linear transformation
  beta = 0.5;      % 
  rrseed = 17; 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f19(x, DIM, ntrial)
% F8F2 sum of Griewank-Rosenbrock 2-D blocks
% last change: 08/12/17
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 19; 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    scale = max(1, sqrt(DIM) / 8.); 
    linearTF = scale * computeRotation(rseed, DIM); 
    Xopt = linearTF' * 0.5*ones(DIM,1) / scale^2;
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = linearTF * x + 0.5;    % rotate/scale

  %----- COMPUTATION core -----
  F2 = 100 * (x(1:end-1,:).^2 - x(2:end,:)).^2 + (1 - x(1:end-1,:)).^2;
  Ftrue = 10 + 10 * sum(F2/4000 - cos(F2), 1) / (DIM - 1); 
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

end % function

%%%%%%%%%%%%%%%%%% weak global structure %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f20(x, DIM, ntrial)
% Schwefel with tridiagonal variable transformation 
% last change: 08/12/17

  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo arrSigns rseed

  funcID = 20; 
  condition = 10;  % for linear transformation 
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = []; % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f21(x, DIM, ntrial)
% Gallagher with 101 Gaussian peaks, condition up to 1000, one global rotation
% CAVE: the difficulty might significantly depend on the
%       random realization for the location of optima. 
% last change: 09/01/26
  persistent Fopt Xopt scales linearTF rotation Xlocal peakvalues
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 21; 
  maxcondition = 1000;  % for linear transformation
  fitvalues = [1.1, 9.1];    % range for raw values, 10 is optimal
  nhighpeaks = 101;   
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    [ignore, idx] = sort(unif(nhighpeaks-1, rseed));  % random permutation
    arrCondition = [sqrt(maxcondition) arrCondition(idx)];
    arrScales = []; 
    for i = 1:nhighpeaks  % generation of cov matrices 
      s = arrCondition(i).^linspace(-0.5, 0.5, DIM); 
      [ignore, idx] = sort(unif(DIM, rseed+1e3*(i-1))); % permutation instead of rotation
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
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = []; % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f22(x, DIM, ntrial)
% Gallagher with 21 Gaussian peaks, condition up to 1000, one global rotation
% CAVE: the difficulty might significantly depend on the
%       random realization for the location of optima. 
% last change: 09/01/26
  persistent Fopt Xopt scales linearTF rotation Xlocal peakvalues
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 22; 
  maxcondition = 1000;  % for linear transformation
  fitvalues = [1.1, 9.1];    % range for raw values, 10 is optimal
  nhighpeaks = 21;   
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    [ignore, idx] = sort(unif(nhighpeaks-1, rseed));  % random permutation
    arrCondition = [maxcondition arrCondition(idx)];
    arrScales = []; 
    for i = 1:nhighpeaks  % generation of cov matrices 
      s = arrCondition(i).^linspace(-0.5, 0.5, DIM); 
      [ignore, idx] = sort(unif(DIM, rseed+1e3*(i-1))); % permutation instead of rotation
      arrScales(i,:) = s(idx); % this is inverse Cov
    end
    % compute peak values, 10 is global optimum
    peakvalues = [10 linspace(fitvalues(1), fitvalues(2), nhighpeaks-1)]; 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    Xlocal = rotation * reshape(9.8*unif(DIM*nhighpeaks, rseed)-4.9, ...
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
      Ftrue{1} = linearTF; 
      Ftrue{2} = rotation; 
    else  % if strcmpi(strinput, 'info')
      Ftrue = []; % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f23(x, DIM, ntrial)
% Katsuura function
% last change: 09/01/26
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 23; 
  condition = 100;  % for linear transformation
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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
    scales = (99.1234/5)^0 * sqrt(condition).^linspace(0, 1, DIM)'; 
    linearTF = diag(scales) * computeRotation(rseed, DIM); 
    % decouple scaling from function definition
    linearTF = rotation * linearTF; % or computeRotation(rseed+1e3, DIM)
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
  Fpen = 1 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
                    % x = rotation * x;  % no scaling here, because it would go to the arrExpo
  % x = monotoneTFosc(x); 
  % idx = x > 0;
  % x(idx) = x(idx).^(1 + arrExpo(idx) .* sqrt(x(idx)));  % smooth in zero
  % x = rotation' * x;   % back-rotation
  % x = arrScales .* x;  % scaling, Xopt should be 0 
  x = linearTF * x;    % rotate/scale

  %----- COMPUTATION core -----
  arr2k = 2.^(1:32);  % TODO: d=32 has an influence? 
  % TODO: can this loop be replaced? 
  for k = 1:POPSI
    arr = x(:,k) * arr2k;  % DIMxd array
    Ftrue(k) = -10/DIM^2 + ...
        10/DIM^2 * prod(1 + (1:DIM)' .* (abs(arr-round(arr)) * arr2k'.^-1)).^(10/DIM^1.2); 
  end
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

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f24(x, DIM, ntrial)
% Lunacek bi-Rastrigin, condition 100
% in PPSN 2008, Rastrigin part rotated and scaled 
% last change: 08/12/17
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 24; 
  condition = 100;  % for linear transformation
  
  %----- CHECK INPUT -----
  if ischar(x) % return Fopt Xopt or linearTF on string argument
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
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

end % function

%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%

function x_opt = computeXopt(seed, DIM)
   % rounded by for digits, but never to zero
   x_opt = 8 * floor(1e4*unif(DIM,seed)')/1e4 - 4;
   idx = (x_opt == 0);
   x_opt(idx) = -1e-5;
end

function B = computeRotation(seed, DIM)
% computes an orthogonal basis
  B = reshape(gauss(DIM*DIM,seed), DIM, DIM);
  for i = 1:DIM
      for j = 1:i-1
	B(:,i) = B(:,i) - B(:,i)'*B(:,j) * B(:,j);
      end
      B(:,i) = B(:,i) / sqrt(sum(B(:,i).^2));
    end
end

function g = monotoneTFosc(f)
% maps [-inf,inf] to [-inf,inf] with different constants
% for positive and negative part
   a = 0.1;
   g = f; 
   idx = (f > 0);
   g(idx) = log(f(idx))/a;
   g(idx) = exp(g(idx) + 0.49*(sin(g(idx)) + sin(0.79*g(idx)))).^a;
   idx = (f < 0);
   g(idx) = log(-f(idx))/a;
   g(idx) = -exp(g(idx) + 0.49*(sin(0.55*g(idx)) + sin(0.31*g(idx)))).^a;
end

%---------- pseudo random number generator ------------
function g = gauss(N, seed)
% gauss(N, seed)
% samples N standard normally distributed numbers
% being the same for a given seed
  r = unif(2*N, seed); % in principle we need only half
  g = sqrt(-2*log(r(1:N))) .* cos(2*pi*r(N+1:2*N));
  if any(g == 0)
    g(g == 0) = 1e-99;
  end
end

function r = unif(N, inseed)
% unif(N, seed)
%    generates N uniform numbers with starting seed

  % initialization
  inseed = abs(inseed);
  if inseed < 1
    inseed = 1;
  end
  aktseed = inseed;
  for i = 39:-1:0
    tmp = floor(aktseed/127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if aktseed < 0
      aktseed = aktseed + 2147483647;
    end
    if i < 32
      rgrand(i+1) = aktseed;
    end
  end
  aktrand = rgrand(1);

  % sample numbers
  r = zeros(1,N); % makes the function ten times faster(!)
  for i = 1:N
    tmp = floor(aktseed/127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if aktseed < 0
      aktseed = aktseed + 2147483647;
    end
    tmp = floor(aktrand / 67108865);
    aktrand = rgrand(tmp+1);
    rgrand(tmp+1) = aktseed;
    r(i) = aktrand/2.147483647e9;
  end
  if any(r == 0)
    warning('zero sampled(?), set to 1e-99');
    r(r == 0) = 1e-99;
  end
end

% for testing and comparing to other implementations, 
%   rand and randn are used only for sampling the noise 
%   should be renamed (ie removed from execution) in the end
%   in 4-D internal rand makes functions 30% faster. 
function res = myrandn(N,M)
   persistent seed
   if isempty(seed)
     seed = 1;
     disp('non-matlab randn');
   else
     seed = seed + 1;  
     if seed > 1e9
       seed = 1;
     end
   end
   res = reshape(gauss(N*M, seed), N, M);
end

function res = myrand(N,M)
   persistent seed
   if isempty(seed)
     seed = 1;
     disp('non-matlab rand');
   else
     seed = seed + 1;
     if seed > 1e9
       seed = 1;
     end
   end
   res = reshape(unif(N*M, seed), N, M); 
end

% qqq
%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = template(x, DIM, ntrial)
% last change: 09/01/29
  persistent Fopt Xopt scales linearTF rotation rseed
  persistent lastSize arrXopt arrScales arrExpo 

  funcID = INPUTTHIS; 
  condition = INPUTTHIS;  % for linear transformation
%  alpha = 1;      % 
%  beta = 0.25;       % 
  rrseed = funcID; 
  
  %----- CHECK INPUT -----
  if ischar(x) % initialize or return Fopt Xopt or linearTF 
    flginputischar = 1;
    strinput = x;
    if nargin < 2 || isempty(DIM)
      DIM = 2;
    end
    x = ones(DIM,1);  % for setting all persistent variables
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
    rseed = rrseed;  % like for ntrial==0
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
    % arrExpo = repmat(beta * linspace(0, 1, DIM)', 1, POPSI); 
    % arrScales = repmat(scales, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  % Fpen(Fpen > 0) = NaN;  % rejection boundary handling 
  % x(xoutside ~= 0) = x(xoutside ~= 0) - xoutside(xoutside ~= 0);  % into [-5, 5]
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  % x = rotation * x;  % no scaling here, because it would go to the arrExpo
  % x = monotoneTFosc(x); 
  % idx = x > 0;
  % x(idx) = x(idx).^(1 + arrExpo(idx) .* sqrt(x(idx)));  % smooth in zero
  % x = rotation' * x;   % back-rotation
  % x = arrScales .* x;  % scaling, Xopt should be 0 
  x = linearTF * x;    % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = 0; 
  Fval = Ftrue;  % without noise

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 
  Fval = FCauchy(Ftrue, 1); 

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

end % function



