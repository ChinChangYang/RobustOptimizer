function estimate_cec11_time
% [lb, ub] = getlimit_messenger;
[lb, ub] = getlimit_cassini2;

startTime = tic;
for t = 1 : 1e5
% Elapsed time is 4.252734 seconds.
% 	cec11_f1(rand(6,1));
% Elapsed time is 15.180237 seconds.
% 	cec11_f7(rand(20,1));
% cec11_f10(0.2 + 0.5 * rand(12, 1));
% Elapsed time is 543.438464 seconds.
% cec11_f12(lb' + rand * (ub' - lb'));
% Elapsed time is 280.334972 seconds.
cec11_f13(lb' + rand * (ub' - lb'));
% Elapsed time is 251.666702 seconds.
end
toc(startTime);
end

