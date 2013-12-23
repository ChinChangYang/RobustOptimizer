function f = twoplane(x, y)
% Two-plane function
% Global min-max value: f(lb,lb) = ?
f = min(3 - 2/10*sum(x) + 3/10*sum(y), ...
	3 + 2/10*sum(x) - 1/10*sum(y));
end
