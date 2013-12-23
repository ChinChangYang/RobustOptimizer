function [x, f, counteval, countsuccess] = de_select(x, u, f, counteval, fitfun, varargin)
% DE_SELECT Selection of Differential Evolution

[~, uNP] = size(u);
[~, xNP] = size(x);
countsuccess = 0;

for i = 1 : uNP
	fui = feval(fitfun, u(:, i));
	counteval = counteval + 1;
	
	if i <= xNP
		if fui < f(i)
			f(i) = fui;
			x(:, i) = u(:, i);
			countsuccess = countsuccess + 1;
		end
	end
end

end