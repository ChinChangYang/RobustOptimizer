function ret = mf5(x)
% MF5 Multiobjective Function No.5
ret = [sphere(x); ellipsoidal(x); differentpowersrot(x); rastrigin(x); ...
	griewankrosenbrock(x)];
end

