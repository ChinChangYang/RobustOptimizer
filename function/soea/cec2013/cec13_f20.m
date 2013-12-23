function f = cec13_f20(x)
% CEC13_F20 CEC'2013 function 20
% Expanded Scaffer's F6 Function
x = reshape(x, numel(x), 1);
f = cec13_func(x, 20) - 600;
if f < 1e-8
	f = 0;
end
end

