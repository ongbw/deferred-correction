function y = vode2_exact(t)

    y = [1./(1+t.^2); -2*t./(1+t.^2).^2];

end