function dfdy = vode2_dfdy(t,y)

    dfdy    = zeros(length(y));
    dfdy(1,2) = 1;
    dfdy(2,1) = -4*y(1)*(1-6*t^2*y(1));

end