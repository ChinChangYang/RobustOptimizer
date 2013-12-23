function f = sainz_f2(x, y)
%SAINZ_F2 Function value of Problem 2 (Sainz et al., 2008)

if nargin == 0
	f = [4.14292612399524,	4.14292612399524; ...
		4.80704852482444,	6.90709922718507];
	return;
end

f = sum(x.^2 + y.^2 + 2 .* x .* y - 20 * x - 20 * y + 100);
end
