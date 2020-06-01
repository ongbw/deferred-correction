function pp = splinederiv(tvals, yvals)

pp = spline(tvals, yvals);

pp.coefs = pp.coefs(:,1:end-1);

for k = 1:(pp.order-1)
    pp.coefs(:,k) = pp.coefs(:,k)*(pp.order-k);   
end
pp.order = pp.order -1;

