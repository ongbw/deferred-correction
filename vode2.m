% this is the right-hand side of problem 1 from Zadunaisky 1976

function dydt = vode2(t,y)

    dydt    = zeros(length(y),1);
    dydt(1) = y(2);
    dydt(2) = -2*y(1)^2*(1-4*t^2*y(1));

end