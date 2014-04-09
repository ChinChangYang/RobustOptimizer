function f = cec14_f10(x)
% CEC14_F10 CEC'2014 function 10
% Shifted Schwefel¡¦s Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 10) - 1000;
end
