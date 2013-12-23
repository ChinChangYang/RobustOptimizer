function r = unif(N, inseed)
% unif(N, seed)
%    generates N uniform numbers with starting seed

  % initialization
  inseed = abs(inseed);
  if inseed < 1
    inseed = 1;
  end
  aktseed = inseed;
  for i = 39:-1:0
    tmp = floor(aktseed/127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if aktseed < 0
      aktseed = aktseed + 2147483647;
    end
    if i < 32
      rgrand(i+1) = aktseed;
    end
  end
  aktrand = rgrand(1);

  % sample numbers
  r = zeros(1,N); % makes the function ten times faster(!)
  for i = 1:N
    tmp = floor(aktseed/127773);
    aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
    if aktseed < 0
      aktseed = aktseed + 2147483647;
    end
    tmp = floor(aktrand / 67108865);
    aktrand = rgrand(tmp+1);
    rgrand(tmp+1) = aktseed;
    r(i) = aktrand/2.147483647e9;
  end
  if any(r == 0)
    warning('zero sampled(?), set to 1e-99');
    r(r == 0) = 1e-99;
  end
end
