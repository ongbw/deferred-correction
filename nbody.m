function ydot = nbody(t, y)

[m,n] = size(y);
nplus = m/4;
nminus = m/4;

%nplus = 200;
%nminus = 200;

ntotal = nplus + nminus;
qplus = 1/nplus;
qminus = -1/nminus;
mplus = 1000/nplus;
mminus = 1/nminus;
d = 0.05;

xplus = y(1:nplus, :);
vplus = y(nplus + 1:2*nplus, :);
xminus = y(2*nplus + 1:2*nplus + nminus, :);
vminus = y(2*nplus + nminus + 1:end, :);

xplusdot = vplus;
xminusdot = vminus;

vplusdot = zeros(size(xminus));
vminusdot = zeros(size(vminus));

for i = 1:nplus
    vplusdot(i, :) = (qplus/mplus)*(sum(qplus*(repmat(xplus(i, :), nplus, 1) - xplus) ./ sqrt((repmat(xplus(i, :), nplus, 1) - xplus).^2 + d^2)) + sum(qminus*(repmat(xplus(i, :), nminus, 1) - xminus) ./ sqrt((repmat(xplus(i, :), nminus, 1) - xminus).^2 + d^2)));
end

for i = 1:nminus
    vminusdot(i, :) = (qminus/mplus)*(sum(qplus*(repmat(xminus(i, :), nplus, 1) - xplus) ./ sqrt((repmat(xminus(i, :), nplus, 1) - xplus).^2 + d^2)) + sum(qminus*(repmat(xminus(i, :), nminus, 1) - xminus) ./ sqrt((repmat(xminus(i, :), nminus, 1) - xminus).^2 + d^2)));
end

ydot = [xplusdot; vplusdot; xminusdot; vminusdot];
