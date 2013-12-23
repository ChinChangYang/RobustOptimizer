function plotfunction(fitfun, L, R, U)
%PLOTFUNCTION Plot a function within L and U. The interval between two
%samples is defined as R.
X = L:R:U;
f = numel(X);
for i = 1 : numel(X)
	f(i) = feval(fitfun, X(i));
end
plot(X, f);
end

