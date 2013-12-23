function plotthreeminmaxfunc
fitfun = 'prtb_quad';
lb = [0.7;-2];
ub = [1.3;4];
X = lb(2):0.01:ub(2);
Y1 = zeros(1, numel(X));
Y2 = Y1;
Y3 = Y1;
for i = 1 : numel(X)
	Y1(i) = feval(fitfun, lb(1), X(i)) + 30;
	Y2(i) = feval(fitfun, (lb(1) + ub(1))/2, X(i)) + 30;
	Y3(i) = feval(fitfun, ub(1), X(i)) + 30;
end
plot(X, Y1, 'k');
hold on;
plot(X, Y2, 'k--');
hold on;
plot(X, Y3, 'k:');
hold off;
legend(sprintf('x = %.2f', lb(1)), ...
	sprintf('x = %.2f', (lb(1) + ub(1))/2), ...
	sprintf('x = %.2f', ub(1)));
xlabel('y');
ylabel('f');
end