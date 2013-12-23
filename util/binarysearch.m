function mid = binarysearch(x, A, left, right)
% Binary Search Algorithm
while left <= right
	mid = round((left + right) / 2);
	if A(mid) < x
		left = mid + 1;
	elseif A(mid) > x
		right = mid - 1;
	else
		return;
	end
end
mid = round((left + right) / 2);
end
