function out = finishoutput(out, X, f, counteval, countiter, varargin)
%FINISHOUTPUT Finish output info
if isempty(out.recordFEs)	
	if ~isempty(varargin)
		for i = 1 : 2 : numel(varargin)
			out.(varargin{i}) = varargin{i + 1};
		end
	end
	
	out.fes = counteval;
	out.G = countiter;
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
			if length(varargin{i + 1}) == 1
				data(iRecordFEs:RecordPoint) = varargin{i + 1};
			else				
				data(:, iRecordFEs:RecordPoint) = ...
					repmat(varargin{i + 1}, 1, numel(iRecordFEs:RecordPoint));
			end
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
out.fmin(iRecordFEs:end) = fmin;
out.fmean(iRecordFEs:end) = mean(f);
out.fstd(iRecordFEs:end) = std(f);
nRemaining = RecordPoint - iRecordFEs + 1;
out.xmin(:, iRecordFEs:RecordPoint) = repmat(xmin, 1, nRemaining);
out.xmean(:, iRecordFEs:RecordPoint) = repmat(mean(X, 2), 1, nRemaining);
out.xstd(:, iRecordFEs:RecordPoint) = repmat(std(X, 0, 2), 1, nRemaining);
out.fes(iRecordFEs:RecordPoint) = out.recordFEs(iRecordFEs:RecordPoint);
out.distancemean(iRecordFEs:RecordPoint) = mean(centroiddistance(X));
out.cond(iRecordFEs:RecordPoint) = condX;
out.angle(iRecordFEs:RecordPoint) = angle;
out.G(iRecordFEs:RecordPoint) = out.recordG(iRecordFEs:RecordPoint);
out.final.X = X;
out.final.f = f;
end
