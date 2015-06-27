function [xmin, fmin, out] = despa(fitfun, lb, ub, maxfunevals, options)
% DESPA A Differential Evolution Algorithm with Success-based Parameter
% Adaptation for CEC2015 Learning-based Optimization
% DESPA(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DESPA(..., options) minimize the function with solver options.
%
% Paper:
% Noor Awad, Mostafa Z. Ali and Robert G. Reynolds, "A Differential
% Evolution Algorithm with Success-based Parameter Adaptation for CEC2015
% Learning-based Optimization," in CEC2015 competition on learning-based
% real-parameter single objective optimization, Sendai International
% Centre, Sendai, Japan, 25-28 May 2015
%
% Special Notes:
% This code is written by Chin-Chang Yang who rewritten the program in
% Matlab. No guarantee of identical behavior to the original algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare Input Arguments                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <= 4
    options = [];
end

if nargin <= 3
    maxfunevals = [];
end

if nargin <= 1
    lb = [];
    ub = [];
end

if nargin <= 0
    fitfun = [];
end

if isempty(fitfun)
    fitfun = 'sphere';
end

if isempty(lb) || isempty(ub)
    lb = -100 * ones(10, 1);
    ub = 100 * ones(10, 1);
end

if isempty(maxfunevals)
    maxfunevals = numel(lb) * 1e4;
end

defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.EarlyStop = 'none';
defaultOptions.usefunevals = inf;
defaultOptions.NP = 4;
defaultOptions.F = 0.5;
defaultOptions.CR = 0.5;
defaultOptions.NPmin = 4;
defaultOptions.NPmax = 450;
defaultOptions.M = 10;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.NP = [];
defaultOptions.initial.MF = [];
defaultOptions.initial.MCR = [];
defaultOptions.initial.iM = [];

options = setdefoptions(options, defaultOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare initial parameters                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common
isDisplayIter   = strcmp(options.Display, 'iter');
RecordPoint     = max(0, floor(options.RecordPoint));
ftarget         = options.ftarget;
usefunevals     = options.usefunevals;

if ~isempty(strfind(options.EarlyStop, 'fitness'))
    EarlyStopOnFitness  = true;
    AutoEarlyStop       = false;
elseif ~isempty(strfind(options.EarlyStop, 'auto'))
    EarlyStopOnFitness  = false;
    AutoEarlyStop       = true;
else
    EarlyStopOnFitness  = false;
    AutoEarlyStop       = false;
end

% Algorithmic constant
NPmin = options.NPmin;
NPmax = options.NPmax;
M = options.M;

% Algorithmic state variables
if ~isempty(options.initial)
    options.initial = setdefoptions(options.initial, defaultOptions.initial);
    
    X   = options.initial.X;
    fx  = options.initial.f;
    NP  = options.initial.NP;
    MF  = options.initial.MF;
    MCR = options.initial.MCR;
    iM	= options.initial.iM;
else
    X = [];
    fx = [];
    NP = [];
    MF = [];
    MCR = [];
    iM	= [];
end

if isempty(NP)
    NP = options.NP;
end

if isempty(MF)
    MF = options.F * ones(1, M);
end

if isempty(MCR)
    MCR = options.CR * ones(1, M);
end

if isempty(iM)
    iM = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Algorithm                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables
D = numel(lb);
counteval = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals, ...
    'NP');

% Initialize contour data
if isDisplayIter
    [XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
end

% Initialize population
if isempty(X)
    X = zeros(D, NPmax);
    for i = 1 : NPmax
        X(:, i) = lb + (ub - lb) .* rand(D, 1);
    end
end

% Evaluation
if isempty(fx)
    fx = nan(1, NPmax);
    for i = 1 : NP
        fx(i) = feval(fitfun, X(:, i));
        counteval = counteval + 1;
    end
end

% Sort
[fx(1 : NP), fidx] = sort(fx(1 : NP));
X(:, 1 : NP) = X(:, fidx);

% Initialize variables
V = X;
U = X;
fu = zeros(1, NPmax);
df = zeros(1, NPmax);
SF = zeros(1, NPmax);
SCR = zeros(1, NPmax);
Chy = cauchyrnd(0, 0.1, NPmax + 10);
iChy = 1;
    
% Algorithmic parameters
r1 = zeros(1, NPmax);
r2 = zeros(1, NPmax);

% Display
if isDisplayIter
    displayitermessages(...
        X(:, 1 : NP), U(:, 1 : NP), fx(1 : NP), countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X(:, 1 : NP), fx(1 : NP), counteval, countiter, ...
    'NP', NP);

% Iteration counter
countiter = countiter + 1;

while true
    % Termination conditions
    outofmaxfunevals = counteval > maxfunevals - NP;
    outofusefunevals = counteval > usefunevals - NP;
    if ~EarlyStopOnFitness && ~AutoEarlyStop
        if outofmaxfunevals || outofusefunevals
            break;
        end
    elseif AutoEarlyStop
        reachftarget = min(fx) <= ftarget;
        TolX = 10 * eps(mean(X(:)));
        solutionconvergence = std(X(:)) <= TolX;
        TolFun = 10 * eps(mean(fx));
        functionvalueconvergence = std(fx(:)) <= TolFun;
        stagnation = countStagnation >= 100;
        
        if outofmaxfunevals || ...
                outofusefunevals || ...
                reachftarget || ...
                solutionconvergence || ...
                functionvalueconvergence || ...
                stagnation
            break;
        end
    elseif EarlyStopOnFitness
        reachftarget = min(fx) <= ftarget;
        
        if outofmaxfunevals || ...
                outofusefunevals || ...
                reachftarget
            break;
        end
    end
    
    % Generate p
    p = floor(1 + max(2, round(NP * 0.11)) * rand(1, NP));
    
    % Generate k
    k = floor(1 + M * rand(1, NP));
    
    % Generate CR
    CR = MCR(k) + 0.1 * randn(1, NP);
    CR(CR < 0) = 0;
    CR(CR > 1) = 1;
    
    % Generate F
    F = zeros(1, NP);
    for i = 1 : NP
        while F(i) <= 0
            F(i) = MF(k(i)) + Chy(iChy);
            iChy = mod(iChy, numel(Chy)) + 1;
        end
    end
    F(F > 1) = 1;
    
    for i = 1 : NP
        % Generate r1
        r1(i) = floor(1 + NP * rand);
        
        % Generate r2
        r2(i) = floor(1 + NP * rand);
        while r1(i) == r2(i)
            r2(i) = floor(1 + NP * rand);
        end
    end
    
    % Mutation
    for i = 1 : NP
        V(:, i) = X(:, i) + F(i) .* (X(:, p(i)) - X(:, i)) + ...
            F(i) .* (X(:, r1(i)) - X(:, r2(i)));
    end
    
    % Crossover
    for i = 1 : NP
        jrand = floor(1 + D * rand);
        for j = 1 : D
            if rand < CR(i) || j == jrand
                U(j, i) = V(j, i);
            else
                U(j, i) = X(j, i);
            end
        end
    end
    
    % Correction for outside of boundaries
    for i = 1 : NP
        for j = 1 : D
            if U(j, i) < lb(j)
                U(j, i) = 0.5 * (lb(j) + X(j, i));
            elseif U(j, i) > ub(j)
                U(j, i) = 0.5 * (ub(j) + X(j, i));
            end
        end
    end
    
    % Display
    if isDisplayIter
        displayitermessages(...
            X, U, fx, countiter, XX, YY, ZZ);
    end
    
    % Evaluation
    for i = 1 : NP
        fu(i) = feval(fitfun, U(:, i));
        counteval = counteval + 1;
    end
    
    % Selection
    FailedIteration = true;
    nS = 0;
    for i = 1 : NP
        if fu(i) < fx(i)
            nS      = nS + 1;
            SF(nS)  = F(i);
            SCR(nS) = CR(i);
            df(nS)  = abs(fu(i) - fx(i));
            
            X(:, i)		= U(:, i);
            fx(i)		= fu(i);
            FailedIteration = false;
        elseif fu(i) == fx(i)
            X(:, i)		= U(:, i);
            fx(i)		= fu(i);
        end
    end
    
    % Update parameter
    if nS > 0
        w = df(1 : nS) ./ sum(df(1 : nS));
        
        MCR(iM) = sum(w .* SCR(1 : nS) .* SCR(1 : nS)) / (sum(w .* SCR(1 : nS)) + eps);
        MF(iM) = sum(w .* SF(1 : nS) .* SF(1 : nS)) / (sum(w .* SF(1 : nS)) + eps);
        
        iM = 1 + mod(iM, M);
    end
    
    % Sort
    [fx(1 : NP), fidx] = sort(fx(1 : NP));
    X(:, 1 : NP) = X(:, fidx);
    
    % Update NP
    if counteval < 0.75 * maxfunevals
        NP = min(NPmax, NP + round((NP + NPmax) * counteval / maxfunevals + NPmin));        
        inan = find(isnan(fx(1 : NP)));
        
        for i = inan
            fx(i) = feval(fitfun, X(:, i));
            counteval = counteval + 1;
        end
    else
        NP = min(NPmax, round((NPmin - NPmax) * counteval / maxfunevals + NPmax));
        
        X = X(:, 1 : NP);
        fx = fx(1 : NP); 
    end
    
    % Sort
    [fx(1 : NP), fidx] = sort(fx(1 : NP));
    X(:, 1 : NP) = X(:, fidx);
    
    % Record
    out = updateoutput(out, X(:, 1 : NP), fx(1 : NP), counteval, countiter, ...
        'NP', NP);
    
    % Iteration counter
    countiter = countiter + 1;
    
    % Stagnation iteration
    if FailedIteration
        countStagnation = countStagnation + 1;
    else
        countStagnation = 0;
    end
end

fmin = fx(1);
xmin = X(:, 1);

final.NP = NP;
final.MF = MF;
final.MCR = MCR;
final.iM = iM;

out = finishoutput(out, X, fx, counteval, countiter, ...
	'final', final, ...
    'NP', NP);
end
