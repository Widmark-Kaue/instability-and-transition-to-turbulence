function [V,lambda, L, F, y, D, D4] = orrSommerfeld(N,alpha, Re, baseFlow)
% Solve Orr-Sommerfeld eigenvalue problem for a specific base flow 

[D, D2, D4, y] = dirichletNeumannBCs(N);
    
% Constants
II = eye(size(D));
alpha2 =  alpha^2;
alpha4 = alpha^4;

% Base Flow
U = baseFlow(y);

% set eigenvalue problem
L = diag(U)*(D2 - alpha2*II) - diag(D2*U) - 1/(1i*alpha*Re) * (D4 - 2*alpha2 * D2 + alpha4*II) ;
F =  D2 - alpha2*II;

% solve eigenvalue problem
[V,lambda]=eig(L,F);
lambda = diag(lambda);
end