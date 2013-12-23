function f = sainz_f1(x, y)
%SAINZ_F1 Function value of Problem 1 (Sainz et al., 2008)

if nargin == 0
	f = [-0.437081694561628,	-0.437081694561628; ...
		-2.55383263412616,		-3.14];
	return;
end

f = sum((cos(y) + cos(2 * y + x)).^2);
end
