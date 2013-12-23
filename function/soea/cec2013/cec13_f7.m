function f = cec13_f7(x)
% CEC13_F7 CEC'2013 function 7
% Rotated Schaffers F7 Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 7) + 800;
if f < 1e-8
	f = 0;
end
end

