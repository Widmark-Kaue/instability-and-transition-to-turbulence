clear all; close all;

set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',16);

figure;

a = 0.1;
b = 0.2;

x0 = 0.3;

    
tspan = 0:0.1:300;
options = odeset('RelTol',1e-4);

[t,q] = ode45(@(t,X) pitchfork(t,X,a,b),tspan,x0,options);

figure;
plot(t,q);
xlabel('t');ylabel('x');
