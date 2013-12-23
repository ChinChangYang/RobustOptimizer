function options = setdefoptions(options, defaultOptions) 
%SETOPTIONS Set options with default options
optionNames = fieldnames(defaultOptions);

for i = 1 : numel(optionNames)
	if ~isfield(options, optionNames{i})
		options.(optionNames{i}) = defaultOptions.(optionNames{i});
	end
end
end

