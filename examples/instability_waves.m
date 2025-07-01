clear all; close all;

%% plot parameters
% set(0,'DefaultLineLinewidth',2);
% set(0,'DefaultAxesFontSize',12);

%% wave parameters

alpha = 1.-0.*1i;  % wavenumber
omega = 0.5+0.1*1i; % frequency

% 1i is the imaginary unit

A = 1; % wave amplitude at t=0

%% spatial grid
lambda = 2*pi/real(alpha);
dx = lambda/64;

x = 0:dx:5*lambda;
nx = length(x);

%% loop in time
T = 2*pi/real(omega); % "period" of oscillation
dt = T/64;

t = 0:dt:5*T;
nt = length(t);

figure;

for i=1:nt
    wave = real(A*exp(1i*(alpha*x-omega*t(i))));
    plot(x,wave);
    xlim([0 max(x)]);ylim([-10*A 10*A]);
    title(['t = ' num2str(t(i))]);
    xlabel('x');ylabel('v(x,y_0,t)');
    pause(0.05);drawnow;   
end
