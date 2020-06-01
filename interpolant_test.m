% testing interpolant

clear

for kk =1:10
    Nx = 40*2^(kk-1);
    x = linspace(0,1,Nx)';
    y = sin(2*pi*x);
    % y = [0; y; 0]; % for clamped spline: gives bad error???

    inty = (cos(2*pi*(x(2:end))) - cos(2*pi*(x(1:end-1))))/(2*pi);

    pp = splineint(x,y);

    quadint = zeros(Nx-1,1);
    for k = 2:(Nx-1)
        quadint(k) = polyval(pp.coefs(k-1,:),x(k)) -  ...
            polyval(pp.coefs(k-1,:),x(k-1));
    end

    err = sum(abs(inty-quadint))/Nx;
    fprintf('Nx = %d, error = %g\n',Nx,err);
end