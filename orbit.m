function [sol] = orbit(t,y)

  mu = 0.012277471;  mup = 1 - mu;
  D1 = (  (y(1) + mu)^2  + y(2).^2 )^(3/2);
  D2 = (  (y(1) - mup)^2 + y(2).^2 )^(3/2);
  sol = zeros(size(y));
  sol(1) = y(3);
  sol(2) = y(4);
  sol(3) = y(1) + 2*y(4) - mup * (y(1) + mu) / D1 ...
           - mu * (y(1) - mup)/D2;
  sol(4) = y(2) - 2*y(3) - mup*y(2) / D1 - mu*y(2)/D2;
  
  %y1_0'= 0.994, y2_0' = 0, y2(0) = 0, y2'(0) = -2.001585
  % tend = 17.065
