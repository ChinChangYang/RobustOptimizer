function g = gauss(N, seed)
% gauss(N, seed)
% samples N standard normally distributed numbers
% being the same for a given seed
  r = unif(2*N, seed); % in principle we need only half
  g = sqrt(-2*log(r(1:N))) .* cos(2*pi*r(N+1:2*N));
  if any(g == 0)
    g(g == 0) = 1e-99;
  end
end
