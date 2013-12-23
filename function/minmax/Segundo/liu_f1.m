function f = liu_f1(x, y)
%LIU_F1 Function value of Problem 4 (Liu et al., 1998)

if nargin == 0
	f = [7.07106781186548,	7.07106781186548; ...
		7.07106781186548,	7.07106781186548;
		7.07106781186548,	0;
		0,					7.07106781186548];
	return;
end

f = -(x(1) + y(1)) .* (x(2) + y(2)) / (1 + x(1) * y(1) + x(2) * y(2));
end
