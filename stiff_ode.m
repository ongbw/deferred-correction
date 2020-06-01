function ydot = stiff_ode(t, y)
epsilon = 1e-5;
ydot = [y(2, :); (1 - y(1, :).^2).*y(2, :)./epsilon - y(1, :)];
