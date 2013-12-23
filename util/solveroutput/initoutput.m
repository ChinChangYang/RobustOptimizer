function out = initoutput(RecordPoint, D, NP, maxfunevals, varargin)
% INITOUTPUT Initialize output info
out.iRecordFEs = 1;
out.recordFEs = linspace(0, 1, RecordPoint) * maxfunevals;
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

if ~isempty(varargin)
	nvarargin = numel(varargin);
	for i = 1 : nvarargin
		out.(varargin{i}) = zeros(1, RecordPoint);
	end
end
end
