clear all; close all; clc;

refL = 1; %% reference length, 1 for delta, 2 for delta*

%set(0,'DefaultLineLinewidth',2);
%set(0,'DefaultAxesFontSize',16);

global beta

%% grid
N = 400; % degree of highest Chebyshev polynomial






%% Orr-Sommerfeld considering viscosity beta = 0.1

[D,y]=cheb(N); %from Trefethen, Spectral Methods in Matlab

H = 20; y = (y+1)*H; D = D/H;

% mean flow
% Blasius boundary layer
delta = 1.;
fact = 1./delta * 5.;

beta = 0.1; % Falkner-Skan; 

options = optimset('TolX',1e-10);
ypp = fminsearch(@findfpp,0.4,options);

f0 = [0 0 ypp];
[xout,yout]=ode45(@blasius,y(N+1:-1:1),f0);

U = yout(N+1:-1:1,2);

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

alpha = 1;
Re = 4500;
II = eye(size(D2));
alpha2 = alpha^2*II;
alpha4 = alpha^4*II;

% define L(v)
L = diag(U)*(D2 - alpha2) - diag(D2*U) - diag(1/(alpha*Re*1i)) * (D4 - 2*alpha2*D2 + alpha4); 
% define F(v)
F = D2 - alpha2;


% solve eigenvalue problem
[V,lambda]=eig(L,F);
lambda = diag(lambda);

%plot real and imaginary parts
figModes = figure;
subplot(1,2,1);
plot(real(lambda),imag(lambda),'bs');
ylim([-2,1]);
title('Eigenspectrum, \alpha = 1 Re = 4500 \beta = 0.1');
yline(0);
xlabel('c_r');
ylabel('c_i');




%% Orr-Sommerfeld considering viscosity beta = -0.1

[D,y]=cheb(N); %from Trefethen, Spectral Methods in Matlab

H = 20; y = (y+1)*H; D = D/H;

% mean flow
% Blasius boundary layer
delta = 1.;
fact = 1./delta * 5.;

beta = -0.1; % Falkner-Skan; 

options = optimset('TolX',1e-10);
ypp = fminsearch(@findfpp,0.4,options);

f0 = [0 0 ypp];
[xout,yout]=ode45(@blasius,y(N+1:-1:1),f0);

U = yout(N+1:-1:1,2);

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

alpha = 1;
Re = 4500;
II = eye(size(D2));
alpha2 = alpha^2*II;
alpha4 = alpha^4*II;

% define L(v)
L = diag(U)*(D2 - alpha2) - diag(D2*U) - diag(1/(alpha*Re*1i)) * (D4 - 2*alpha2*D2 + alpha4); 
% define F(v)
F = D2 - alpha2;


% solve eigenvalue problem
[V,lambda]=eig(L,F);
lambda = diag(lambda);

%plot real and imaginary parts
figure(figModes);
subplot(1,2,2);
plot(real(lambda),imag(lambda),'bs');
ylim([-2,1]);
title('Eigenspectrum, \alpha = 1 Re = 4500 \beta = -0.1');
yline(0);
xlabel('c_r');
ylabel('c_i');








%% Neutral curve and critical Reynolds number

alphaVec = .5:0.05:1.2;
ReVec = [.5:.05:10]*1000 ;
alphaZero = zeros(length(alphaVec),0);
ReZero = zeros(length(alphaVec),0);

for i=1:length(alphaVec)
    alpha = alphaVec(i);
    for j=1:length(ReVec)
        Re = ReVec(j);
        II = eye(size(D2));
        alpha2 = alpha^2*II;
        alpha4 = alpha^4*II;
        
        % define L(v)
        L = diag(U)*(D2 - alpha2) - diag(D2*U) - diag(1/(alpha*Re*1i)) * (D4 - 2*alpha2*D2 + alpha4); 
        % define F(v)
        F = D2 - alpha2;
        
        
        % solve eigenvalue problem
        [V,lambda]=eig(L,F);
        lambda = diag(lambda);

        pos = find(imag(lambda) > 0.01); %row position of unstable modes in matrix

        if isempty(pos)
            alphaZero(i) = alphaVec(i);
            ReZero(i) = ReVec(j);
        end
    end
end

plot(ReZero, alphaZero);
xline(min(ReZero));
xlim([1000 10000]);
ylim([0 1.2])
title('Neutral Curve (\omega_i=0)');
xlabel('Re');
ylabel('\alpha');
txt = strcat(num2str(min(ReZero)), '\rightarrow');
text(min(ReZero),0.6,txt,'HorizontalAlignment','right');






