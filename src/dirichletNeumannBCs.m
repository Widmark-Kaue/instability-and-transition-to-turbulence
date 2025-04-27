function [D,D2,D4,y] = dirichlet_neumann_BCs(N)
% Compute Chebyshev differentiation matrix with Dirichlet and Neumann
% boundary conditions
% 
[D, y] = cheb(N);  %from Trefethen, Spectral Methods in Matlab
D2 = D*D; 
S = diag([0; 1./(1-y(2:N).^2); 0]);
D4 = (diag(1-y.^2)*D^4 - 8*diag(y)*D^3 - 12*D^2)*S;

y = y(2:N);
D = D(2:N,2:N);
D2 = D2(2:N,2:N);
D4 = D4(2:N,2:N);
end