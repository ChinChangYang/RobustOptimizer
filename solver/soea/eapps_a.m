function [xmin, fmin, out] = eapps_a(fitfun, lb, ub, maxfunevals, options)
% EAPPS_A Evolutionary Algorithm with Probability-based Parent-Selecting
% Framework (Variant A)
% EAPPS_A(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% EAPPS_A(..., options) minimize the function with solver options.

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
defaultOptions.NP = 150;
defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];

options = setdefoptions(options, defaultOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare initial parameters                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F               = options.F;

isDisplayIter   = strcmp(options.Display, 'iter');
RecordPoint     = max(0, floor(options.RecordPoint));
ftarget         = options.ftarget;

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

if ~isempty(options.initial)
    options.initial = setdefoptions(options.initial, defaultOptions.initial);
    X = options.initial.X;
    fx = options.initial.f;
else
    X = [];
    fx = [];
end

if isempty(X)
    NP = options.NP;
else
    [~, NP] = size(X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Algorithm                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables
D = numel(lb);
counteval = 0;
countiter = 1;
countStagnation = 0;
out = initoutput(RecordPoint, D, NP, maxfunevals);

% Initialize contour data
if isDisplayIter
    [XX, YY, ZZ] = advcontourdata(D, lb, ub, fitfun);
end

% Initialize population
if isempty(X)
    X = zeros(D, NP);
    for i = 1 : NP
        X(:, i) = lb + (ub - lb) .* rand(D, 1);
    end
end

% Evaluation
if isempty(fx)
    fx = zeros(1, NP);
    for i = 1 : NP
        fx(i) = feval(fitfun, X(:, i));
        counteval = counteval + 1;
    end
end

% Sort
[fx, fidx] = sort(fx);
X = X(:, fidx);

% Initialize variables
V = X;
U = X;
fu = zeros(1, NP);
r1 = zeros(1, NP);
r2 = zeros(1, NP);
r3 = zeros(1, NP);
CR = zeros(1, NP);
iCR = zeros(1, NP);
qCR = 1 / 11 * ones(1, 11);

% Display
if isDisplayIter
    displayitermessages(...
        X, U, fx, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, fx, counteval, countiter);

% Iteration counter
countiter = countiter + 1;

while true
    % Termination conditions
    outofmaxfunevals = counteval > maxfunevals - NP;
    if ~EarlyStopOnFitness && ~AutoEarlyStop
        if outofmaxfunevals
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
                reachftarget || ...
                solutionconvergence || ...
                functionvalueconvergence || ...
                stagnation
            break;
        end
    elseif EarlyStopOnFitness
        reachftarget = min(fx) <= ftarget;
        
        if outofmaxfunevals || ...
                reachftarget
            break;
        end
    end
    
    for i = 1 : NP
        % Generate r1
        r1(i) = floor(1 + NP * rand);
        
        % Generate r2
        r2(i) = floor(1 + NP * rand);
        while r1(i) == r2(i)
            r2(i) = floor(1 + NP * rand);
        end
        
        % Generate r3
        r3(i) = floor(1 + NP * rand);
        while r1(i) == r3(i) || r2(i) == r3(i)
            r3(i) = floor(1 + NP * rand);
        end
        
        % Generate CR
        iCR(i) = sum(rand >= cumsum([0, qCR]));
        CR(i) = (iCR(i) - 1) / 10;
    end
    
    % Mutation
    for i = 1 : NP
        V(:, i) = X(:, r1(i)) + F .* (X(:, r2(i)) - X(:, r3(i)));
    end
    
    % Binominal Crossover
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
    for i = 1 : NP
        if fu(i) < fx(i)
            X(:, i)		= U(:, i);
            fx(i)		= fu(i);
            
            qCR(iCR(i))   = qCR(iCR(i)) * 2;
            FailedIteration = false;
        else
            qCR(iCR(i))   = qCR(iCR(i)) * 0.9170;
        end
    end
    
    % Sort
    [fx, fidx] = sort(fx);
    X = X(:, fidx);
    
    % Update parameter probability
    qCR = qCR./sum(qCR);
    
    % Record
    out = updateoutput(out, X, fx, counteval, countiter);
    
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

out = finishoutput(out, X, fx, counteval, countiter);
end
