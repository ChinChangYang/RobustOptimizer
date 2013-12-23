function f = functionsampler(x, options)
%FUNCTIONSAMPLER Function sampler
if nargin <= 1
	options = [];
end
defaultOptions.fitfun = 'cloplog';
defaultOptions.T = 100;
options = setdefoptions(options, defaultOptions);
fitfun = options.fitfun;
T = options.T;
S = zeros(1, T);
for i = 1 : T
	S(i) = feval(fitfun, x);
end
f = mean(S);
end
