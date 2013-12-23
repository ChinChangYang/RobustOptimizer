function f = cec11_f12(x)
persistent initial_flag MGADSMproblem %#ok<USENS>
if isempty(initial_flag)
	load messengerfull.mat;
	initial_flag = 1;
end
f = messengerfull(x, MGADSMproblem);
end

