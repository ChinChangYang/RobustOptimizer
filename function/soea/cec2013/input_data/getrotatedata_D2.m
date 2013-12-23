fileID = fopen('M_D2.txt');
tline = fgetl(fileID);
A = [];
while ischar(tline)
	A = [A, sscanf(tline, '%e')];
	tline = fgetl(fileID);
end
fclose(fileID);