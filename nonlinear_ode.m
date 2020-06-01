function ydot = nonlinear_ode(t, y)
ydot = [-y(2, :) + y(1,:).*(1- y(1, :).^2 - y(2,:).^2);
        y(1, :) + +3*y(2,:).*(1- y(1, :).^2 - y(2,:).^2)];

