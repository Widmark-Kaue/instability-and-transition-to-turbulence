function bl=blasius(t,X)
global beta
f = X(1); fp = X(2); fpp = X(3);
bl=zeros(size(X));
bl(1) = fp; bl(2) = fpp; bl(3) = -0.5*(f*fpp + beta*(1-fp^2)) ;