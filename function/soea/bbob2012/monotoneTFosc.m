function g = monotoneTFosc(f)
% maps [-inf,inf] to [-inf,inf] with different constants
% for positive and negative part
   a = 0.1;
   g = f; 
   idx = (f > 0);
   g(idx) = log(f(idx))/a;
   g(idx) = exp(g(idx) + 0.49*(sin(g(idx)) + sin(0.79*g(idx)))).^a;
   idx = (f < 0);
   g(idx) = log(-f(idx))/a;
   g(idx) = -exp(g(idx) + 0.49*(sin(0.55*g(idx)) + sin(0.31*g(idx)))).^a;
end
