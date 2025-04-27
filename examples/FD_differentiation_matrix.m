clear all; close all;

%% grid
N = 40; % number of gridpoints
x0 = -1; xf = 1; Lx = xf-x0; dx = Lx/(N-1);
x = (x0:dx:xf).';

%% second-order centered differences for first derivative

um = ones(N,1);

D = full(spdiags([-1/2*um 0*um 1/2*um],-1:1,N,N));
D(1,:)=0;D(1,1) = -1; D(1,2)=1;
D(N,:)=0;D(N,N-1) = -1; D(N,N)=1; % first-order differences for boundary points
D = D/dx;

%% test derivative
f = exp(x.^2);
df = D*f;
dfan = 2*x.*(exp(x.^2));
plot(x,df,'bs',x,dfan,'k-');
xlabel('x');ylabel('f''');
legend('Numerical derivative','Analytical derivative');

max(abs(df-dfan)) %maximum error

