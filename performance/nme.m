function ret = nme(me)
% NME Normalized mean error
[nA, nf] = size(me);
if nA == 1 || nf == 1
	if any(me < 0)
		% Exception: negative error detected
		% Resolution: shift errors to be larger than zero
		me = me - 2 * min(me);
	end
	ret = (me + eps) / (min(me) + eps);
	return;
end

ret = me;
for i = 1 : nf
	if any(me(:, i) < 0)
		% Exception: negative error detected
		% Resolution: shift errors to be larger than zero
		me(:, i) = me(:, i) - 2 * min(me(:, i));
	end
	ret(:, i) = (me(:, i) + eps) ./ (min(me(:, i)) + eps);
end
end
