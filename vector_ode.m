function ydot = vector_ode(t, y)
ydot = [t.*y(2, :) + y(1, :); -t.*y(1, :) + y(2, :)];
