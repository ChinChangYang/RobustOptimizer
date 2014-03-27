load('errors_cec2013_201403261747.mat')
A1			= reshape(results, 28, 51)';
A1_MEAN		= mean(results, 3)';
A1_STD		= std(results, [], 3)';

load('errors_cec2013_201403271144.mat')
A2			= reshape(results, 28, 51)';
A2_MEAN		= mean(results, 3)';
A2_STD		= std(results, [], 3)';

w			= ranksumtest(A1, A2);
POSITIVE	= sum(w=='+');
EQUAL		= sum(w=='=');
NEGATIVE	= sum(w=='-');
