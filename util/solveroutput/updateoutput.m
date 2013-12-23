function out = updateoutput(out, X, f, counteval, varargin)
%UPDATEOUTPUT Update output info
if isempty(out.recordFEs)
	return;
end

if counteval >= out.recordFEs(out.iRecordFEs)
	[fmin, fminidx] = min(f);
	xmin = X(:, fminidx);
	distance = sampledistance(X);
	C = cov(X');
	[B, ~] = eig(C);
	condX = cond(C);
	angle = acos(B(1));
	if angle >= 0.5 * pi
		angle = angle - 0.5 * pi;
	end
	while counteval >= out.recordFEs(out.iRecordFEs)
		i = out.iRecordFEs;
		out.fmin(i) = fmin;
		out.fmean(i) = mean(f);
		out.fstd(i) = std(f);
		out.xmin(:, i) = xmin;
		out.xmean(:, i) = mean(X, 2);
		out.xstd(:, i) = std(X, 0, 2);
		out.fes(i) = counteval;
		out.distancemin(i) = min(distance);
		out.distancemax(i) = max(distance);
		out.distancemean(i) = mean(distance);
		out.distancemedian(i) = median(distance);
		out.cond(i) = condX;
		out.angle(i) = angle;
		out.iRecordFEs = out.iRecordFEs + 1;
		
		if ~isempty(varargin)
			for j = 1 : 2 : numel(varargin)
				data = out.(varargin{j});
				data(i) = varargin{j + 1};
				out.(varargin{j}) = data;
			end
		end
		
		if out.iRecordFEs > numel(out.fmin)
			break;
		end
	end
end
end
