function [Fval, Ftrue] = bbob12_f9(x)
% Rosenbrock, rotated
% last change: 08/12/10
  persistent Fopt Xopt linearTF 
  persistent lastSize rseed

  funcID = 9; 
  
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
    scale = max(1, sqrt(DIM) / 8.); 
    linearTF = scale * computeRotation(rseed, DIM); 
    Xopt = linearTF' * 0.5*ones(DIM,1) / scale^2;
  end
  % DIM- and POPSI-dependent initializations of DIMxPOPSI matrices
  if isempty(lastSize) || lastSize.DIM ~= DIM || lastSize.POPSI ~= POPSI
    lastSize.POPSI = POPSI; 
    lastSize.DIM = DIM; 
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
    else  % if strcmpi(strinput, 'info')
      Ftrue = []; % benchmarkinfos(funcID); 
      Fval = Fopt;
    end
  end
end
