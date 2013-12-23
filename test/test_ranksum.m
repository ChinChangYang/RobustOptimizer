function test_ranksum
T = 1000;
N = 30;
a = 0.1;
h = zeros(1, T);

for t = 1 : T
	[~, h(t), ~] = ranksum(randn(1, N), randn(1, N), a);
end

disp(mean(h));
end

