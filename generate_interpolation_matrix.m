function L = generate_interpolation_matrix(x,xp)
% function L = generate_interpolation_matrix(x,xp)
%
% This function generates the interpolation matrix for a vector of
% quadrature nodes x specified on [0,1]. for various xp's in [0,1]/ 
%
% Specifically, L(j,:) gives the interpolation weights for
% approximating f(xp(j))

m = length(x);
mp = length(xp);

L = zeros(mp,m);

for j=1:mp
  for k=1:m
    coeff = 1;
    for l=1:k-1
      % form interpolating polynomial
      coeff = conv(coeff,[1,-x(l)])/(x(k)-x(l));
    end
    for l=k+1:m
      % form interpolating polynomial
      coeff = conv(coeff,[1,-x(l)])/(x(k)-x(l));
    end
    L(j,k) = polyval(coeff,xp(j));
  end
end
