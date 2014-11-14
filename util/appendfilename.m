function newmetafilename = appendfilename(metafilename, filename)
load(metafilename);
filenames{end+1} = filename; %#ok<NASGU>
metafiledate = datestr(now, 'yyyymmddHHMM');
newmetafilename = sprintf('filenames_%s.mat', metafiledate);
save(newmetafilename, 'filenames');
end

