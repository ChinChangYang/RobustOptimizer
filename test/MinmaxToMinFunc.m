function f = MinmaxToMinFunc( x )
%MinmaxToMinFunc Transform the constraint minimax problems to min problems
%in order to find the true optima.
%
% Note: Find the minimum of this function by:
% > [xmin,fmin,out]=jadebin('MinmaxToMinFunc',lb,ub,1e6)

% P1-Sainz-2008: lb=-0.5; ub=-0.4
f = abs(sainz_f1(x, -3.14) - sainz_f1(x, x*(x+6.28)));

% P2-Sainz-2008: lb=4; ub=4.5
% f = abs(sainz_f2(x, 3 + sqrt(4 - (x - 5)^2)) - ...
% 	sainz_f2(x, 3 + sqrt(16 - (x - 5)^2)));

% P3-Lu-2008: lb=-1; ub=0
% f = lu_f1(x, 5);
end

