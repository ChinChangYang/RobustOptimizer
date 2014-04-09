function f = cec14_f16(x)
% CEC14_F16 CEC'2014 function 16
% Shifted and Rotated Expanded Scaffer¡¦s F6 Function
x = reshape(x, numel(x), 1);
f = cec14_func(x, 16) - 1600;
end
