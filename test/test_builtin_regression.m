function test_builtin_regression
%TEST_BUILTIN_REGRESSION Test built-in regression in MATLAB
x = -1:0.0001:1;
y = zeros(1, numel(x));
for i = 1 : numel(x)
	y(i) = cloplog(x(i));
end
p = polyfit(x,y,2);
f = polyval(p,x);
plot(x,y,'x',x,f,'-');
axis([-1 1 0 1]);
end
