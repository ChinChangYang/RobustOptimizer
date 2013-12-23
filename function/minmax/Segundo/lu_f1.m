function f = lu_f1(x, y)
%LU_F1 Function value of Problem 3 (LU et al., 2008)

if nargin == 0
	f = [-0.887343319848318; ...
		5];
	return;
end

f = sum(sin(x).^2 - x .* cos(y) + 2 * sin(x) - cos(y).^2 + y - 1);
end
