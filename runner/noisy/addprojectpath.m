function addprojectpath
%ADDPROJECTPATH Add project pathes

% Goto project root
cd ..
cd ..

% Get path information
currentPath = cd;
functionPath = sprintf('%s/function', currentPath);
performancePath = sprintf('%s/performance', currentPath);
solverPath = sprintf('%s/solver', currentPath);
utilPath = sprintf('%s/util', currentPath);
runnerPath = sprintf('%s/runner', currentPath);

% Set path
addpath(genpath(functionPath));
addpath(genpath(performancePath));
addpath(genpath(solverPath));
addpath(genpath(utilPath));
addpath(genpath(runnerPath));
end

