function options = setDEoptions(x, t)
options.NP = round(x(1, t)) + round(x(13, t));
options.F = x(2, t);
options.CR = x(3, t);
options.ER = x(4, t);
options.p = x(5, t);
options.H = round(x(6, t));
options.Q = round(x(7, t));
options.Ar = round(x(8, t));
options.cw = x(9, t);
options.erw = x(10, t);
options.CRmin = x(11, t);
options.CRmax = x(11, t) + x(12, t);
options.NPmin = round(x(13, t));
options.crw = x(14, t);
options.fw = x(15, t);
options.EarlyStop = 'auto';
options.RecordFEsFactor = inf;
end

