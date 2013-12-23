function [X, Y, Z, C, h] = drawcontour( fitfun, L, R, U, X, Y, Z )
if nargin < 1
	error('fitfun is not defined.');
end

if nargin < 4
	L = -5 * ones(1, 2);
	U = 5 * ones(1, 2);
	R = 0.05 * ones(1, 2);
end

if nargin < 7
	X = [];
	Y = [];
	Z = [];
end

if isscalar(L)
    L = L * ones(1, 2);
end

if isscalar(R)
    R = R * ones(1, 2);
end

if isscalar(U)
    U = U * ones(1, 2);
end

if isempty(X) || isempty(Y) || isempty(Z)
	X = zeros(length(L(1):R(1):U(1)), length(L(2):R(2):U(2)));
	Y = zeros(size(X));
	Z = zeros(size(X));
	Xidx = 1;
	for i = L(1):R(1):U(1)
		Yidx = 1;
		for j = L(2):R(2):U(2)
			X(Xidx, Yidx) = i;
			Y(Xidx, Yidx) = j;
			Z(Xidx, Yidx) = feval(fitfun, [i; j]);
			Yidx = Yidx + 1;
		end
		Xidx = Xidx + 1;
	end
	Z = log(Z - min(min(Z))+1);
end

[C, h] = contour(X, Y, Z, 30, 'LineWidth', 2);
xlabel('x1');
ylabel('x2');
end
