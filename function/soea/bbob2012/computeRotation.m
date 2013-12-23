function B = computeRotation(seed, DIM)
% computes an orthogonal basis
  B = reshape(gauss(DIM*DIM,seed), DIM, DIM);
  for i = 1:DIM
      for j = 1:i-1
	B(:,i) = B(:,i) - B(:,i)'*B(:,j) * B(:,j);
      end
      B(:,i) = B(:,i) / sqrt(sum(B(:,i).^2));
    end
end
