function f = cec14_f18(x)
% CEC14_F18 CEC'2014 function 18
% Hybrid Function 2
% N = 3
% p = [0.3, 0.3, 0.4]
% g1 : Bent Cigar Function f2
% g2 : HGBat Function f12
% g3: Rastrigin¡¦s Function f8
x = reshape(x, numel(x), 1);
f = cec14_func(x, 18) - 1800;
end
