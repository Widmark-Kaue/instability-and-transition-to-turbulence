close all; clear all;

%% grid
N = 120; % number of gridpoints
y0 = -1; yf = 1; Ly = yf-y0; dy = Ly/(N-1);
y = (y0:dy:yf).';

%% second-order centered differences for second derivative

um = ones(N,1);

D2 = full(spdiags([1*um -2*um 1*um],-1:1,N,N));
D2(1,:)=D2(2,:);D2(N,:)=D2(N-1,:); % first-order differences for boundary points
D2 = D2/(dy^2);

%% set eigenvalue problem
L = D2;
F = -eye(size(L));

%% set boundary conditions v=0 at both ends
L(1,:) = 0; L(1,1) = 1; F(1,:) = 0;
L(N,:) = 0; L(N,N) = 1; F(N,:) = 0;

%% solve eigenvalue problem
[V,lambda]=eig(L,F);
lambda = diag(lambda);

%% sort eigenvalues
lambda=sort(lambda);

%% analytical eigenvalues
n = (1:length(lambda)).';
lambdaAn = n.^2*pi^2/4;

error = abs(lambda-lambdaAn);

normerror = error/abs(lambdaAn); % error normalised by eigenvalue magnitude

figure;
semilogy(n,error,'k.-');
line([0 max(n)],[0.01 0.01],'LineStyle','--'); % draw line marking an error of 0.01
xlabel('n');
ylabel('Absolute error of nth eigenvalue');

%% look at eigenfunctions
j = 10; %% look at jth eigenfunction
if(mod(j,2)==0)
    Van = sin(n(j)*pi/2*y);
else
    Van = cos(n(j)*pi/2*y);
end
figure;
plot(y,V(:,j),'k.',y,Van,'k-');
xlabel('y'); title(['Eigenfunction number ' int2str(j)]);
legend('Numerical','Analytical');
