function metafilename3 = mergefilenames(metafilename1, metafilename2)
load(metafilename2);
temp = filenames; %#ok<NODEF>
n = numel(filenames);
load(metafilename1);
for i = 1 : n
	filenames{end+1} = temp{i};  %#ok<AGROW>
end
metafilename3 = sprintf('filenames_%s', datestr(now, 'yyyymmddHHMM'));
save(metafilename3, 'filenames');
end
