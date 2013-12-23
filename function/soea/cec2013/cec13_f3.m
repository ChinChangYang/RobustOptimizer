function f = cec13_f3(x)
% CEC13_F3 CEC'2013 function 3
% Rotated Bent Cigar Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 3) + 1200;
if f < 1e-8
	f = 0;
end
end

