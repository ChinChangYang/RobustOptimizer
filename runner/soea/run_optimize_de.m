D = 10;
maxfes = 6000;
for i = 1 : 15
	fnum = i;
	optimize_de(D, fnum, maxfes);
end

% D = 30;
% maxfes = 2000;
% for i = 1 : 15
% 	fnum = i;
% 	optimize_de(D, fnum, maxfes);
% end
% 
% D = 50;
% maxfes = 1200;
% for i = 1 : 15
% 	fnum = i;
% 	optimize_de(D, fnum, maxfes);
% end
% 
% D = 100;
% maxfes = 600;
% for i = 1 : 15
% 	fnum = i;
% 	optimize_de(D, fnum, maxfes);
% end
