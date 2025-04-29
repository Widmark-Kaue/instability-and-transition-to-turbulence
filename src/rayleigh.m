function [V,lambda, L, F, y, D, D2] = rayleigh(N,alpha, baseFlow, H)
%Solve Rayleigh eigenvalue problem for a specific base flow 

% grid and derivatives
[D,y]=cheb(N);
D2 = D*D;

% Stretche domain from [-1:1] to [-H:H]
if nargin < 4
    H = 1;
end
y = y*H;
D = D/H;
D2 = D2/(H^2);

% Constants
alpha2 =  alpha^2*eye(size(D2));
U = baseFlow(y);

% set eigenvalue problem
L = diag(U)*(D2 - alpha2) - diag(D2*U) ;
F = D2 - alpha2;

% set boundary conditions v=0 at both ends
L(1,:) = 0; L(1,1) = 1; F(1,:) = 0;
L(N,:) = 0; L(N,N) = 1; F(N,:) = 0;

% solve eigenvalue problem
[V,lambda]=eig(L, F);
lambda = diag(lambda);

end