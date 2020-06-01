function x = chebyshevnodes(a,b,n)

% Given an interval [a, b] and the number of desired
% points n, compute the Chebyshev nodes on [a, b]

k = [1:n];
x = cos( (2*k - 1)/(2*n)*pi);

x = (a+b)/2 + (b-a)/2*x;
x = sort(x);

 