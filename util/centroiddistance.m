function ret = centroiddistance(x)
%CENTROIDDISTANCE Distance between each of x and the centroid of x
c = mean(x, 2);
[~, NP] = size(x);
ret = zeros(1, NP);
for i = 1 : NP
	ret(i) = norm(x(:, i) - c);
end
end
