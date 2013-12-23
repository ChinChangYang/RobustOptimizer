function x_opt = computeXopt(seed, DIM)
   % rounded by for digits, but never to zero
   x_opt = 8 * floor(1e4*unif(DIM,seed)')/1e4 - 4;
   idx = (x_opt == 0);
   x_opt(idx) = -1e-5;
end
