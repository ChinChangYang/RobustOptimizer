function plot_all_maxminmax_func

if matlabpool('size') == 0
	matlabpool
end

parfor i = 1 : 64
	fitfun = sprintf('maxminmax_f%d', i);
	plot_maxminmax_func(fitfun);
end

matlabpool close
end
