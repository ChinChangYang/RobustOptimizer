function d = sampledistance(x)
%SAMPLEDISTANCE Sample the distances between each pair of vectors
[~, NP] = size(x);
if NP <= 1
	d = 0;
	return;
elseif NP == 2
	d = norm(x(:,1) - x(:,2));
	return;
elseif NP == 3
	d = zeros(NP, 1);
	d(1) = norm(x(:,1) - x(:,2));
	d(2) = norm(x(:,2) - x(:,3));
	d(3) = norm(x(:,1) - x(:,3));
	return;
end

d = zeros(NP, 1);
for i = 1 : NP
	r1 = floor(1 + rand * NP);
	r2 = floor(1 + rand * NP);
	while r1 == r2
		r2 = floor(1 + rand * NP);
	end
	d(i) = norm(x(:,r1) - x(:,r2));
end
end
