function plot_all_minmax_func

if matlabpool('size') == 0
	matlabpool
end

parfor i = 1 : 8
	fitfun = sprintf('fminmax_f%d', i);
	plot_minmax_func(fitfun);
end

matlabpool close
end
