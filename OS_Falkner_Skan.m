clear all; close all;

refL = 1; %% reference length, 1 for delta, 2 for delta*

set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',16);

global beta

%% grid
N = 100; % degree of highest Chebyshev polynomial
[D,y]=cheb(N); %from Trefethen, Spectral Methods in Matlab

H = 20; y = (y+1)*H; D = D/H;

% mean flow
% Blasius boundary layer
delta = 1.;
fact = 1./delta * 5.;

beta = 0.0; % Falkner-Skan; 
% put beta=0 for Blasius
% put beta>0 for favourable pressure gradient
% put beta<0 for adverse pressure gradient

options = optimset('TolX',1e-10);
ypp = fminsearch(@findfpp,0.4,options);

f0 = [0 0 ypp];
[xout,yout]=ode45(@blasius,y(N+1:-1:1),f0);

U = yout(N+1:-1:1,2);

%% find boundary layer thickness and normalise y
if(refL == 1)
delta = interp1(U,y,0.99);
elseif(refL == 2)
delta = -trapz(y,(1-U));    
end
y = y/delta; D=D*delta;
fact = H/delta;


Dold = D; D2old = D*D;
dU = D*U;
ddU = D*D*U;

%% derivative operators
[D,y]=cheb(N);
D2 = D*D; D2 = D2(2:N,2:N);
S = diag([0; 1./(1-y(2:N).^2); 0]);
D4 = (diag(1-y.^2)*D^4 - 8*diag(y)*D^3 - 12*D^2)*S;
D4=D4(2:N,2:N); y = y(2:N);
U = U(2:N);
dU = dU(2:N);
ddU = ddU(2:N);
D = D(2:N,2:N);

y = (y+1)*fact;
D = D/fact;D2=D2/fact^2; D4=D4/fact^4;


figure;
plot(U,y,dU,y,ddU,y);
ylabel('y/\delta'); legend('U','U''','U''''');grid on;
