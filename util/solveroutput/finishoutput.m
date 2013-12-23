function out = finishoutput(out, X, f, counteval, varargin)
%FINISHOUTPUT Finish output info
if isempty(out.recordFEs)	
	if ~isempty(varargin)
		for i = 1 : 2 : numel(varargin)
			out.(varargin{i}) = varargin{i + 1};
		end
	end
	
	out.fes = counteval;
	out.final.X = X;
	out.final.f = f;
	return;
end

iRecordFEs = out.iRecordFEs;
RecordPoint = numel(out.fmin);

if ~isempty(varargin)
	for i = 1 : 2 : numel(varargin)
		if isfield(out, varargin{i})
			data = out.(varargin{i});
			data(iRecordFEs:RecordPoint) = varargin{i + 1};
			out.(varargin{i}) = data;
		else
			out.(varargin{i}) = varargin{i + 1};
		end
	end
end

C = cov(X');
[B, ~] = eig(C);
condX = cond(C);
angle = acos(B(1));
if angle >= 0.5 * pi
	angle = angle - 0.5 * pi;
end
[fmin, fminidx] = min(f);
xmin = X(:, fminidx);
distance = sampledistance(X);
out.fmin(iRecordFEs:end) = fmin;
out.fmean(iRecordFEs:end) = mean(f);
out.fstd(iRecordFEs:end) = std(f);
nRemaining = RecordPoint - iRecordFEs + 1;
out.xmin(:, iRecordFEs:RecordPoint) = repmat(xmin, 1, nRemaining);
out.xmean(:, iRecordFEs:RecordPoint) = repmat(mean(X, 2), 1, nRemaining);
out.xstd(:, iRecordFEs:RecordPoint) = repmat(std(X, 0, 2), 1, nRemaining);
out.fes(iRecordFEs:RecordPoint) = counteval;
out.distancemin(iRecordFEs:RecordPoint) = min(distance);
out.distancemax(iRecordFEs:RecordPoint) = max(distance);
out.distancemean(iRecordFEs:RecordPoint) = mean(distance);
out.distancemedian(iRecordFEs:RecordPoint) = median(distance);
out.cond(iRecordFEs:RecordPoint) = condX;
out.angle(iRecordFEs:RecordPoint) = angle;
out.final.X = X;
out.final.f = f;
end
