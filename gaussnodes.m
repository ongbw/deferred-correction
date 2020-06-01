function x = gaussnodes(a,b,n)

% Given an interval [a, b] and the number of desired
% points n, compute the Legendre-Gauss nodes on [a, b]

j = 1:n-1;
beta = 0.5 ./ sqrt(1 - (2*j).^(-2));
T = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(T);
x = diag(D);
x = sort(x);
x = (b - a)/2 * x' + (a + b)/2;