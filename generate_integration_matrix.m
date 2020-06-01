function S = generate_integration_matrix(x,xp)
% function S = generate_integration_matrix(x,xp)
%
% This function generates the integration matrix for a vector of
% quadrature nodes x specified on [0,1]. for various xp's in [0,1]/
%
% Specifically, S(j,:) gives the quadrature weights for
% approximating int(0,xp(j))

m = length(x);
mp = length(xp);
S = zeros(mp,m);

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
    % integrate interpolating polynomial
    coeff_int = polyint(coeff); 
    S(j,k) = polyval(coeff_int,xp(j)) - ...
             polyval(coeff_int,0);
  end
end
