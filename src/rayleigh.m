function [V,lambda, L, F, y, D, D2] = rayleigh(N,alpha, baseFlow, H, h)
%Solve Rayleigh eigenvalue problem for a specific base flow 

% grid and derivatives
[D,y]=cheb(N-1);
D2 = D*D;

% Stretche domain from [-1:1] to [-H:H]
switch nargin
    case 3
        H = 1;
        h = 0;
    case 4
        h = 0;
end
y = (y+h)*H;
D = D/H;
D2 = D2/(H^2);

% Apply boundary conditions v = 0 at both ends
% y = y(2:N+1);
% D = D(2:N+1, 2:N+1);
% D2 = D2(2:N+1, 2:N+1);


% Constants
alpha2 =  alpha^2*eye(size(D2));
U = baseFlow;
if isa(baseFlow, 'function_handle')
    U = baseFlow(y);
end


    


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