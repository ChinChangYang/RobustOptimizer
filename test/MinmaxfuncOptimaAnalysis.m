lb1 = -5;
lb2 = -5;
ub1 = 5;
ub2 = 5;
fitfun = 'lu_f1';
D=1;
options.nonlcon = 'lu_c1';
[XX, YY, ZZ, CC] = cminmaxcontourdata(D, lb1, ub1, lb2, ub2, fitfun, options);
X=[-0.887343319848318; ...
		5];
U=[-0.887343319848318; ...
		5];
f=0;
countiter=1;
dispconitermsg(X, U, f, countiter, XX, YY, ZZ, CC);
title('Problem P3 (Lu et al., 2008)')