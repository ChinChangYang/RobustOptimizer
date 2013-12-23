% VAL = BENCHMARKSNOISY(X, FUNCID)
% VAL = BENCHMARKSNOISY(X, STRFUNC)
%    Input: 
%       X -- solution column vector or matrix of column vectors
%       FUNCID -- number of function to be executed with X as input
%       STRFUNC -- function as string to be executed with X as input
%    Output: function value(s) of solution(s)
%    Examples: 
%      F = BENCHMARKSNOISY([1 2 3]', 117); 
%      F = BENCHMARKSNOISY([1 2 3]', 'f101'); 
% 
% NBS = BENCHMARKSNOISY() 
% NBS = BENCHMARKSNOISY('FunctionIndices') 
%    Output: 
%      NBS -- array of valid benchmark function numbers, 
%             presumably 101:130
%
% FHS = BENCHMARKSNOISY('handles')
%    Output: 
%      FHS -- cell array of function handles
%
% STR = BENCHMARKSNOISY('info', FUNC_NB)
%    Output: 
%      STR -- function description string of function FUNC_NB
%             FUNC_NB == SSCANF(STR, '%d').  
% Examples:
%    fhs = benchmarks('handles');  
%    f = fh{1}(x);  % evaluates x on the sphere with mod noise, f101
%    f = feval(fh{1}, x); % dito
%    fidx = benchmarks('FunctionIndices');
%
% see also: functions FGENERIC, BENCHMARKS, BENCHMARKINFOS

% INTERFACE OF BENCHMARK FUNCTIONS
% see benchmarks.m 

% Authors (copyright 2009): Nikolaus Hansen, Raymond Ros, Steffen Finck
%    Version = '$Revision: 3162 $'
%    Last Modified: $Date: 2011-02-10 16:24:43 +0100 (Thu, 10 Feb 2011) $


%%%-------------------------------------------------------------%%%
function [res, res2] = benchmarksnoisy(x, strfun, DIM) 
%
  Nfcts = 30;

  % return valid function IDs (ie numbers)
  if nargin < 1 || ( ...
      ischar(x) && (strcmpi(x, 'FunctionIndices') || strcmpi(x, 'noisyFunctionIndices')))
    res = [];
    for i = 101:200
      try  % exist does not work under Octave
        eval(['f' num2str(i) '([1; 2]);']);
        res(end+1) = i;
      catch
        if i < 100 + Nfcts
          disp(sprintf('execution of function %d produces an error', i));
        end
      end
    end 
    if length(res) > Nfcts
      disp(res);
      error([num2str(length(res)) ' > 30 functions f1, f2,...,f24, f101,...' ...
             ' are defined. Try to clear the workspace']);
    elseif length(res) < Nfcts
      disp(res);
      error([num2str(length(res)) ' < 30 functions are defined. ' ...
             'There might be an execution error ']);
    end
  elseif isnumeric(x)
    if isnumeric(strfun)
      res = feval(['f' num2str(strfun)], x);
    else
      res = feval(strfun, x);
    end
  % return function handles
  elseif ischar(x) && strcmpi(x, 'handles')
    res = {};
    for i = 1:100
      try  % exist does not work under Octave
        eval(['f' num2str(i+100) '([1; 2]);']);
        res{i} = str2func(['f' num2str(i+100)]); % eval(['@f' num2str(i)]);
      catch
        if i < Nfcts
          disp(sprintf('execution of function %d produces an error', i+100));
        end
      end
    end
    return;
  elseif ischar(x) && strcmpi(x, 'test')
    error('to be implemented');
    % TODO
    for DIM = [2,3,5,10,20,40]
      % generate test points also outside the domain 
      % pop = generateTestPoints
      for ifun = benchmarksnoisy('FunctionIndices')
        % evaluate all functions at all the points, using fgeneric
      end
    end
    % dump the result 
  elseif ischar(x) && nargin < 2  % TODO: simplify
    res = benchmarkinfos(); 
    res = res(101:end);
    if nargin > 1
      res2 = res;
      res = res2{strfun};
    end
  elseif ischar(x) && strcmpi(x, 'fopt')
    if isnumeric(strfun)
      res = feval(['f' num2str(strfun)], x);
    else
      res = feval(strfun, x);
    end
  elseif ischar(x) 
    res = benchmarkinfos(strfun); 
  else
    error('not a valid case of input');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%% NOISY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f101(x, DIM, ntrial)
% sphere with moderate Gauss noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 101; 
  rrseed = 1; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2, 1);

  %----- NOISE -----
  Fval = FGauss(Ftrue, 0.01); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f102(x, DIM, ntrial)
% sphere with moderate uniform noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 102; 
  rrseed = 1; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2, 1);

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.01 * (0.49 + 1/DIM), 0.01); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f103(x, DIM, ntrial)
% sphere with moderate Cauchy noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 103; 
  rrseed = 1; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2, 1);

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 0.01, 0.05); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f104(x, DIM, ntrial)
% Rosenbrock non-rotated with moderate Gauss noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 104; 
  rrseed = 8; 
  
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
    Xopt =1* 0.75 * computeXopt(rseed, DIM); 
    scales = max(1, sqrt(DIM) / 8.); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = scales * x;  % rotate/scale while assuming that Xopt==0 
  x = x + 1;  % shift zero to factual optimum 1

  %----- COMPUTATION core -----
  Ftrue = 1e2 * sum((x(1:end-1, :).^2 - x(2:end, :)).^2, 1) ...
          + sum((x(1:end-1, :) - 1).^2, 1);

  %----- NOISE -----
  Fval = FGauss(Ftrue, 0.01); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f105(x, DIM, ntrial)
% Rosenbrock non-rotated with moderate uniform noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 105; 
  rrseed = 8; 
  
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
    Xopt =1* 0.75 * computeXopt(rseed, DIM); 
    scales = max(1, sqrt(DIM) / 8.); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = scales * x;  % rotate/scale while assuming that Xopt==0 
  x = x + 1;  % shift zero to factual optimum 1

  %----- COMPUTATION core -----
  Ftrue = 1e2 * sum((x(1:end-1, :).^2 - x(2:end, :)).^2, 1) ...
          + sum((x(1:end-1, :) - 1).^2, 1);

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.01 * (0.49 + 1/DIM), 0.01); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f106(x, DIM, ntrial)
% Rosenbrock non-rotated with moderate Cauchy noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 106; 
  rrseed = 8; 
  
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
    Xopt =1* 0.75 * computeXopt(rseed, DIM); 
    scales = max(1, sqrt(DIM) / 8.); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = scales * x;  % rotate/scale while assuming that Xopt==0 
  x = x + 1;  % shift zero to factual optimum 1

  %----- COMPUTATION core -----
  Ftrue = 1e2 * sum((x(1:end-1, :).^2 - x(2:end, :)).^2, 1) ...
          + sum((x(1:end-1, :) - 1).^2, 1);

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 0.01, 0.05); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function


%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f107(x, DIM, ntrial)
% sphere with  Gauss noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 107; 
  rrseed = 1; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2, 1);

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f108(x, DIM, ntrial)
% sphere with uniform noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 108; 
  rrseed = 1; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2, 1);

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f109(x, DIM, ntrial)
% sphere with Cauchy noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 109; 
  rrseed = 1; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 

  %----- COMPUTATION core -----
  Ftrue = sum(x.^2, 1);

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f110(x, DIM, ntrial)
% Rosenbrock non-rotated with Gauss noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 110; 
  rrseed = 8; 
  
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
    Xopt =1* 0.75 * computeXopt(rseed, DIM); 
    scales = max(1, sqrt(DIM) / 8.); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = scales * x;  % rotate/scale while assuming that Xopt==0 
  x = x + 1;  % shift zero to factual optimum 1

  %----- COMPUTATION core -----
  Ftrue = 1e2 * sum((x(1:end-1, :).^2 - x(2:end, :)).^2, 1) ...
          + sum((x(1:end-1, :) - 1).^2, 1);

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f111(x, DIM, ntrial)
% Rosenbrock non-rotated with uniform noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 111; 
  rrseed = 8; 
  
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
    Xopt =1* 0.75 * computeXopt(rseed, DIM); 
    scales = max(1, sqrt(DIM) / 8.); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = scales * x;  % rotate/scale while assuming that Xopt==0 
  x = x + 1;  % shift zero to factual optimum 1

  %----- COMPUTATION core -----
  Ftrue = 1e2 * sum((x(1:end-1, :).^2 - x(2:end, :)).^2, 1) ...
          + sum((x(1:end-1, :) - 1).^2, 1);

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f112(x, DIM, ntrial)
% Rosenbrock non-rotated with Cauchy noise
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 112; 
  rrseed = 8; 
  
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
    Xopt =1* 0.75 * computeXopt(rseed, DIM); 
    scales = max(1, sqrt(DIM) / 8.); 
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = scales * x;  % rotate/scale while assuming that Xopt==0 
  x = x + 1;  % shift zero to factual optimum 1

  %----- COMPUTATION core -----
  Ftrue = 1e2 * sum((x(1:end-1, :).^2 - x(2:end, :)).^2, 1) ...
          + sum((x(1:end-1, :) - 1).^2, 1);

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f113(x, DIM, ntrial)
% step-ellipsoid with gauss noise, condition 100
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 113; 
  condition = 100; 
  alpha = 10;  % inner rounding procedure
  rrseed = 7; 
  
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
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = linearTF * x;  % rotate while assuming that Xopt==0 
  x1 = x(1,:);
  idx = abs(x) > 0.5;
  x(idx) = round(x(idx));
  x(~idx) = round(alpha*x(~idx))/alpha; 
  x = rotation * x;  % rotate while assuming that Xopt==0 

  %----- COMPUTATION core -----
  Ftrue = 0.1 * max(1e-4*abs(x1), scales' * x.^2);

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f114(x, DIM, ntrial)
% step-ellipsoid with uniform noise, condition 100
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 114; 
  condition = 100; 
  alpha = 10;  % inner rounding procedure
  rrseed = 7; 
  
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
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = linearTF * x;  % rotate while assuming that Xopt==0 
  x1 = x(1,:);
  idx = abs(x) > 0.5;
  x(idx) = round(x(idx));
  x(~idx) = round(alpha*x(~idx))/alpha; 
  x = rotation * x;  % rotate while assuming that Xopt==0 

  %----- COMPUTATION core -----
  Ftrue = 0.1 * max(1e-4*abs(x1), scales' * x.^2);

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f115(x, DIM, ntrial)
% step-ellipsoid with Cauchy noise, condition 100
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 115; 
  condition = 100; 
  alpha = 10;  % inner rounding procedure
  rrseed = 7; 
  
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
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = linearTF * x;  % rotate while assuming that Xopt==0 
  x1 = x(1,:);
  idx = abs(x) > 0.5;
  x(idx) = round(x(idx));
  x(~idx) = round(alpha*x(~idx))/alpha; 
  x = rotation * x;  % rotate while assuming that Xopt==0 

  %----- COMPUTATION core -----
  Ftrue = 0.1 * max(1e-4*abs(x1), scales' * x.^2);

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f116(x, DIM, ntrial)
% ellipsoid with Gauss noise, monotone x-transformation, condition 1e4
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 116; 
  condition = 1e4;  
  rrseed = 10; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;  % no scaling here, because it would go to the arrExpo
  x = monotoneTFosc(x);

  %----- COMPUTATION core -----
  Ftrue = scales' * x.^2; 

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function
%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f117(x, DIM, ntrial)
% ellipsoid with uniform noise, monotone x-transformation, condition 1e4
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 117; 
  condition = 1e4;  
  rrseed = 10; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;  % no scaling here, because it would go to the arrExpo
  x = monotoneTFosc(x);

  %----- COMPUTATION core -----
  Ftrue = scales' * x.^2; 
  Fval = Ftrue;  % without noise

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function
%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f118(x, DIM, ntrial)
% ellipsoid with Cauchy noise, monotone x-transformation, condition 1e4
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed

  funcID = 118; 
  condition = 1e4;  
  rrseed = 10; 
  
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
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
    arrXopt = repmat(Xopt, 1, POPSI); 
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;  % no scaling here, because it would go to the arrExpo
  x = monotoneTFosc(x);

  %----- COMPUTATION core -----
  Ftrue = scales' * x.^2; 
  Fval = Ftrue;  % without noise

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f119(x, DIM, ntrial)%
% sum of different powers with Gauss noise, between x^2 and x^6
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 119; 
  alpha = 4; 
  rrseed = 14; 
  
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
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;    % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = sqrt(sum(abs(x).^arrExpo, 1)); 

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f120(x, DIM, ntrial)%
% sum of different powers with uniform noise, between x^2 and x^6
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 120; 
  alpha = 4; 
  rrseed = 14; 
  
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
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;    % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = sqrt(sum(abs(x).^arrExpo, 1)); 

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f121(x, DIM, ntrial)%
% sum of different powers with seldom Cauchy noise, between x^2 and x^6
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 121; 
  alpha = 4; 
  rrseed = 14; 
  
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
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = x - arrXopt;  % shift optimum to zero 
  x = rotation * x;    % rotate/scale

  %----- COMPUTATION core -----
  Ftrue = sqrt(sum(abs(x).^arrExpo, 1)); 

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f122(x, DIM, ntrial)
% Schaffers F7 with Gauss noise, with asymmetric non-linear transformation, condition 10
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 122; 
  condition = 10;  % for linear transformation
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
  Fpen = 100 * sum(xoutside.^2, 1);  
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

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f123(x, DIM, ntrial)
% Schaffers F7 with uniform noise, asymmetric non-linear transformation, condition 10
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 123; 
  condition = 10;  % for linear transformation
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
  Fpen = 100 * sum(xoutside.^2, 1);  
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

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f124(x, DIM, ntrial)
% Schaffers F7 with seldom Cauchy noise, asymmetric non-linear transformation, condition 10
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 124; 
  condition = 10;  % for linear transformation
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
  Fpen = 100 * sum(xoutside.^2, 1);  
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

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f125(x, DIM, ntrial)
% F8F2 sum of Griewank-Rosenbrock 2-D blocks with Gauss noise 
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 125; 
  rrseed = 19; 
  
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
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = linearTF * x + 0.5;    % rotate/scale

  %----- COMPUTATION core -----
  F2 = 100 * (x(1:end-1,:).^2 - x(2:end,:)).^2 + (1 - x(1:end-1,:)).^2;
  Ftrue = 1 + sum(F2/4000 - cos(F2), 1) / (DIM - 1); 

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f126(x, DIM, ntrial)
% F8F2 sum of Griewank-Rosenbrock 2-D blocks with uniform noise 
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 126; 
  rrseed = 19; 
  
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
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = linearTF * x + 0.5;    % rotate/scale

  %----- COMPUTATION core -----
  F2 = 100 * (x(1:end-1,:).^2 - x(2:end,:)).^2 + (1 - x(1:end-1,:)).^2;
  Ftrue = 1 + sum(F2/4000 - cos(F2), 1) / (DIM - 1); 

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f127(x, DIM, ntrial)
% F8F2 sum of Griewank-Rosenbrock 2-D blocks with seldom Cauchy noise 
% last change: 09/01/04
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 127; 
  rrseed = 19; 
  
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
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = linearTF * x + 0.5;    % rotate/scale

  %----- COMPUTATION core -----
  F2 = 100 * (x(1:end-1,:).^2 - x(2:end,:)).^2 + (1 - x(1:end-1,:)).^2;
  Ftrue = 1 + sum(F2/4000 - cos(F2), 1) / (DIM - 1); 

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function


%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f128(x, DIM, ntrial)
% Gallagher with 101 Gaussian peaks with Gauss noise, condition up to 1000, one global rotation
% last change: 09/01/26
% CAVE: the difficulty might significantly depend on the
%       random realization for the location of optima. 
  persistent Fopt Xopt scales linearTF rotation Xlocal peakvalues
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 128; 
  maxcondition = 1000;  % for linear transformation
  fitvalues = [1.1, 9.1];  % range for raw values, 10 is optimal
  nhighpeaks = 101;  % TODO: use 49 for performance reason? 
  rrseed = 21; 
  
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
    Xlocal = rotation * reshape(10.*unif(DIM*nhighpeaks, rseed)-5., ...
                                DIM, nhighpeaks); 
    % global optimum not too close to boundary
    Xlocal(:,1) = 0.8 * Xlocal(:,1); 
    Xopt = rotation' * Xlocal(:,1);
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = rotation * x;  

  %----- COMPUTATION core -----
  fac = -0.5/DIM; % ??xopt + alpha*ones versus alpha becomes invariant like this 
  f = NaN*ones(nhighpeaks,POPSI); 
  if POPSI < 0.5*nhighpeaks  
    for k = 1:POPSI
      xx = repmat(x(:,k), 1, nhighpeaks) - Xlocal;
      f(:,k) = peakvalues' .* exp(fac*(sum(arrScales'.*xx.^2, 1)))';
    end
  else
    for i = 1:nhighpeaks  % CAVE: POPSI=1e4 gets out of memory
      xx = (x - repmat(Xlocal(:,i), 1, POPSI)); % repmat could be done once
      f(i,:) = peakvalues(i) * exp(fac*(arrScales(i,:)*xx.^2)); 
    end
  end
  % f is in [0,10], 10 is best
  Ftrue = monotoneTFosc(10 - max(f, [], 1)).^2;  

  %----- NOISE -----
  Fval = FGauss(Ftrue, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f129(x, DIM, ntrial)
% Gallagher with 101 Gaussian peaks with uniform noise, condition up to 1000, one global rotation
% last change: 09/01/26
% CAVE: the difficulty might significantly depend on the
%       random realization for the location of optima. 
  persistent Fopt Xopt scales linearTF rotation Xlocal peakvalues
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 129; 
  maxcondition = 1000;  % for linear transformation
  fitvalues = [1.1, 9.1];  % range for raw values, 10 is optimal
  nhighpeaks = 101;  % TODO: use 49 for performance reason? 
  rrseed = 21; 
  
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
    Xlocal = rotation * reshape(10.*unif(DIM*nhighpeaks, rseed)-5., ...
                                DIM, nhighpeaks); 
    % global optimum not too close to boundary
    Xlocal(:,1) = 0.8 * Xlocal(:,1); 
    Xopt = rotation' * Xlocal(:,1);
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = rotation * x;  

  %----- COMPUTATION core -----
  fac = -0.5/DIM; % ??xopt + alpha*ones versus alpha becomes invariant like this 
  f = NaN*ones(nhighpeaks,POPSI); 
  if POPSI < 0.5*nhighpeaks  
    for k = 1:POPSI
      xx = repmat(x(:,k), 1, nhighpeaks) - Xlocal;
      f(:,k) = peakvalues' .* exp(fac*(sum(arrScales'.*xx.^2, 1)))';
    end
  else
    for i = 1:nhighpeaks  % CAVE: POPSI=1e4 gets out of memory
      xx = (x - repmat(Xlocal(:,i), 1, POPSI)); % repmat could be done once
      f(i,:) = peakvalues(i) * exp(fac*(arrScales(i,:)*xx.^2)); 
    end
  end
  % f is in [0,10], 10 is best
  Ftrue = monotoneTFosc(10 - max(f, [], 1)).^2;  

  %----- NOISE -----
  Fval = FUniform(Ftrue, 0.49 + 1/DIM, 1); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = f130(x, DIM, ntrial)
% Gallagher with 101 Gaussian peaks with seldom Cauchy noise, condition up to 1000, one global rotation
% last change: 09/01/26
% CAVE: the difficulty might significantly depend on the
%       random realization for the location of optima. 
  persistent Fopt Xopt scales linearTF rotation Xlocal peakvalues
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = 130; 
  maxcondition = 1000;  % for linear transformation
  fitvalues = [1.1,9.1];  % range for raw values, 10 is optimal
  nhighpeaks = 101;  % TODO: use 49 for performance reason? 
  rrseed = 21; 
  
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
    Xlocal = rotation * reshape(10.*unif(DIM*nhighpeaks, rseed)-5., ...
                                DIM, nhighpeaks); 
    % global optimum not too close to boundary
    Xlocal(:,1) = 0.8 * Xlocal(:,1); 
    Xopt = rotation' * Xlocal(:,1);
  end

  %----- BOUNDARY HANDLING -----
  xoutside = max(0, abs(x) - 5) .* sign(x); 
  Fpen = 100 * sum(xoutside.^2, 1);  % penalty
  Fadd = Fadd + Fpen; 

  %----- TRANSFORMATION IN SEARCH SPACE -----
  x = rotation * x;  

  %----- COMPUTATION core -----
  fac = -0.5/DIM; % ??xopt + alpha*ones versus alpha becomes invariant like this 
  f = NaN*ones(nhighpeaks,POPSI); 
  if POPSI < 0.5*nhighpeaks  
    for k = 1:POPSI
      xx = repmat(x(:,k), 1, nhighpeaks) - Xlocal;
      f(:,k) = peakvalues' .* exp(fac*(sum(arrScales'.*xx.^2, 1)))';
    end
  else
    for i = 1:nhighpeaks  % CAVE: POPSI=1e4 gets out of memory
      xx = (x - repmat(Xlocal(:,i), 1, POPSI)); % repmat could be done once
      f(i,:) = peakvalues(i) * exp(fac*(arrScales(i,:)*xx.^2)); 
    end
  end
  % f is in [0,10], 10 is best
  Ftrue = monotoneTFosc(10 - max(f, [], 1)).^2;  

  %----- NOISE -----
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function

% qqq

%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%

function x_opt = computeXopt(seed, DIM)
   % rounded by for digits, but never to zero
   x_opt = 8 * floor(1e4*unif(DIM,seed)')/1e4 - 4;
   idx = x_opt == 0;
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
    g(g==0) = 1e-99;
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
    r(r==0) = 1e-99;
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

function Fval = FGauss(Ftrue, beta)
  POPSI = size(Ftrue, 2);
  Fval = Ftrue .* exp(beta * randn(1, POPSI)); % with gauss noise
  TOL = 1e-8; 
  Fval = Fval + 1.01*TOL; 
  idx = Ftrue < TOL;
  Fval(idx) = Ftrue(idx); 
end

function Fval = FUniform(Ftrue, alpha, beta)
  % alpha = 0.49 + 1/DIM;  % alpha * rand must always be smaller than one 
  % beta = 1;              % smaller is easier
  POPSI = size(Ftrue, 2);
  Fval = rand(1,POPSI).^beta .* Ftrue ...
         .* max(1, (10^9 ./ (Ftrue+1e-99)).^(alpha * rand(1, POPSI)));
  TOL = 1e-8; 
  Fval = Fval + 1.01*TOL; 
  idx = Ftrue < TOL;
  Fval(idx) = Ftrue(idx); 
end

function Fval = FCauchy(Ftrue, alpha, p)
  % Cauchy with median 1e3*alpha and with p=0.2, zero otherwise
  % P(Cauchy > 1,10,100,1000) = 0.25, 0.032, 0.0032, 0.00032
  % alpha = 1;
  POPSI = size(Ftrue, 2);
  Fval = Ftrue + alpha * max(0, 1e3 + (rand(1, POPSI) < p) .* ... 
                          randn(1, POPSI) ./ (abs(randn(1, POPSI))+1e-199)); 
  TOL = 1e-8; 
  Fval = Fval + 1.01*TOL; 
  idx = Ftrue < TOL;
  Fval(idx) = Ftrue(idx); 
end

% qqq
%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%-------------------------------------------------------------%%%
function [Fval, Ftrue] = template(x, DIM, ntrial)
% last change: 09/01/03
  persistent Fopt Xopt scales linearTF rotation
  persistent lastSize arrXopt arrScales arrExpo rseed 

  funcID = INPUTTHIS; 
  condition = INPUTTHIS;  % for linear transformation
%  alpha = 1;      % 
%  beta = 0.25;       % 
  rrseed = INPUTTHIS;  % original function number
  
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
    Fopt =1* min(1000, max(-1000, round(100*100*gauss(1,rndseed)/gauss(1,rndseed+1))/100));
  end 
  Fadd = Fopt;  % value to be added on the "raw" function value
  % DIM-dependent initialization
  if isempty(lastSize) || lastSize.DIM ~= DIM  
    Xopt =1* computeXopt(rndseed, DIM); % function ID is seed for rotation 
    rotation = computeRotation(rndseed+1e6, DIM); 
    scales = sqrt(condition).^linspace(0, 1, DIM)'; 
    linearTF = diag(scales) * computeRotation(rndseed, DIM); 
    % decouple scaling from function definition
    % linearTF = rotation * linearTF; % or computeRotation(rndseed+1e3, DIM)
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
  Fval = FCauchy(Ftrue, 1, 0.2); 

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
      Ftrue = [];  % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end

end % function



