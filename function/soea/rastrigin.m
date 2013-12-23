function f = rastrigin(x)
% RASTRIGIN Rastrigin function
% Properties Highly multimodal function with a comparatively regular
% structure for the placement of the optima. 
%
% * roughly 1e+D local optima 
%
M = numel(x);
x = reshape(x, 1, M);
f = 10 * (M - sum(cos(2*pi*x))) + sum(x.^2);
end

