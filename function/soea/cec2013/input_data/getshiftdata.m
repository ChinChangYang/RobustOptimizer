fileID = fopen('shift_data.txt');
tline = fgetl(fileID);
A = [];
while ischar(tline)
	A = [A, sscanf(tline, '%e')];
	tline = fgetl(fileID);
end
fclose(fileID);