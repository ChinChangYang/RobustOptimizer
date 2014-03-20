function f = cec14_f24(x)
% CEC14_F24 CEC'2014 function 24
% Composition Function 2
% N = 3
% sigma = [20, 20, 20]
% lambda = [1, 1, 1]
% bias = [0, 100, 200]
% g1: Schwefel's Function F10・
% g2: Rotated Rastrigin・s Function F9・
% g3: Rotated HGBat Function F14・
x = reshape(x, numel(x), 1);
f = cec14_func(x, 24) - 2400;
if f < 1e-8
	f = 0;
end
end
