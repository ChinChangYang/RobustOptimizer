function f = cec11_f13(x)
persistent initial_flag MGADSMproblem %#ok<USENS>
if isempty(initial_flag)
	load cassini2.mat;
	initial_flag = 1;
end
f = cassini2(x, MGADSMproblem);
end

