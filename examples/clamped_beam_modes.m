close all; clear all;

%% grid
N = 180; % degree of highest Chebyshev polynomial
[D,y]=cheb(N); %from Trefethen, Spectral Methods in Matlab
D2 = D*D; D2 = D2(2:N,2:N);
S = diag([0; 1./(1-y(2:N).^2); 0]);
D4 = (diag(1-y.^2)*D^4 - 8*diag(y)*D^3 - 12*D^2)*S;
D4=D4(2:N,2:N); y = y(2:N);
D = D(2:N,2:N);
%% set eigenvalue problem

II = eye(size(D));

L = D4;
F = II;
% eigenvalue problem for free vibration of a clamped-clamped
% Euler-Bernoulli beam

%% solve eigenvalue problem
[V,lambda]=eig(L,F);
lambda = diag(lambda);

%% sort eigenvalues
[aux,pos] = sort(lambda);
lambda = lambda(pos); V = V(:,pos);

%% plot eigenmodes

figure;
for i=1:12
subplot(3,4,i);plot(y,V(:,i));title(['Mode ' int2str(i)]);
end