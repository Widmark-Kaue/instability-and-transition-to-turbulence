close all; clear all;

%% grid
N = 360; % degree of highest Chebyshev polynomial
[D,y]=cheb(N); %from Trefethen, Spectral Methods in Matlab
D2 = D*D; D2 = D2(2:N,2:N);
S = diag([0; 1./(1-y(2:N).^2); 0]);
D4 = (diag(1-y.^2)*D^4 - 8*diag(y)*D^3 - 12*D^2)*S;
D4=D4(2:N,2:N); y = y(2:N);
D = D(2:N,2:N);


%% Strech Domain
% H = 20;
% y = y*H;
% D2 = D2/(H^2); D4 = D4/(H^4); 

%% Constants
II = eye(size(D));
alpha = 1;
alpha2 =  alpha^2;
alpha4 = alpha^4;
Re = 10000; 
U = 1 - y.^2;

%% Plot U velocity field
figure
plot(U, y, 'k');
% ylim([-2,2])
xlabel('U')
ylabel('y')
%% set eigenvalue problem
L = diag(U)*(D2 - alpha2*II) - diag(D2*U) - 1/(1i*alpha*Re) * (D4 - 2*alpha2 * D2 + alpha4*II) ;
F = D2 - alpha2*II;

% L = D4;
% F = II;
% eigenvalue problem for free vibration of a clamped-clamped
% Euler-Bernoulli beam

%% solve eigenvalue problem
[V,lambda]=eig(L,F);
lambda = diag(lambda);

%% sort eigenvalues
[aux,pos] = sort(lambda);
lambda = lambda(pos); V = V(:,pos);


%% Plot eingspectrum
figure
plot(real(lambda), imag(lambda), 'ko')
ylim([-2, 1])
xlabel('c_r')
ylabel('c_i')
grid()
%% plot eigenmodes

figure;
for i=1:12
subplot(3,4,i);plot(y,V(:,i));title(['Mode ' int2str(i)]);
end