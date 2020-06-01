function pp = splineint(tvals, yvals)

pp = spline(tvals, yvals);
[m,n] = size(pp.coefs);

for k = 1:pp.order
    pp.coefs(:,k) = pp.coefs(:,k)/(pp.order+1-k);
end
pp.order = pp.order + 1;
pp.coefs(:,pp.order) = zeros(m,1);

