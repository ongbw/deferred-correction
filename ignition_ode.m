function ydot = ignition_ode(t, y)
ydot = y.^2 - y.^3;
