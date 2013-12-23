function ret = mf2( x )
% MF2 Multiobjective Function No.2
ret = [differentpowersrot(x); griewankrosenbrock(x)];
end

