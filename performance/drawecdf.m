function drawecdf(me)
% DRAWECDF Draw empirical cumulative distribution function of mean errors
normalizedErr = nme(me);
[nA, nf] = size(normalizedErr);
ecdf_x = zeros(nA, nf);
ecdf_y = (1 : nf) / nf;
for i = 1 : nA
	ecdf_x(i, :) = sort(normalizedErr(i, :));
	linespec = getlinespec(i);
	semilogx(ecdf_x(i, :), ecdf_y, linespec);
	hold on;
end
xlabel('NME');
ylabel('ECDF');
end
