function simulate_issue76
A = [1,2,3,4,5,6,7,8,9,10];
for i = 1 : 10
	r = 10 * rand;
	index = binarysearch(r, A, 1, 10);
	fprintf('r = %f\n', r);
	fprintf('index = %f\n', index);
end
end
