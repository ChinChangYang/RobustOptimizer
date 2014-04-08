function out = initoutput(RecordPoint, D, NP, maxfunevals, varargin)
% INITOUTPUT Initialize output info
out.iRecordFEs = 1;
out.recordFEs = round(linspace(0, 1, RecordPoint) * maxfunevals);
out.recordG = round(linspace(0, 1, RecordPoint) * maxfunevals / NP);
out.fmin = inf(1, RecordPoint);
out.fmean = inf(1, RecordPoint);
out.fstd = inf(1, RecordPoint);
out.xmean = inf(D, RecordPoint);
out.xmin = inf(D, RecordPoint);
out.xstd = inf(D, RecordPoint);
out.fes = zeros(1, RecordPoint);
out.distance = zeros(NP, RecordPoint);
out.cond = zeros(1, RecordPoint);
out.angle = zeros(1, RecordPoint);
out.bestever.fmin = Inf;
out.G = zeros(1, RecordPoint);

if ~isempty(varargin)
	nvarargin = numel(varargin);
	for i = 1 : nvarargin
		if isequal(varargin{i}, 'FC') || ...
				isequal(varargin{i}, 'MF') || ...
				isequal(varargin{i}, 'MCR')
			out.(varargin{i}) = zeros(NP, RecordPoint);
		else
			out.(varargin{i}) = zeros(1, RecordPoint);
		end
	end
end
end
