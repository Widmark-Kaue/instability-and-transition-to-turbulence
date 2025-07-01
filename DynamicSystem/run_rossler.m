clear all; close all;

set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',16);

figure;

a = 0.1;
b = 0.2;
c = 5.7;

x0 = rand(1,3);
figure;
hold on;

    

tspan = 0:0.1:1000;
options = odeset('RelTol',1e-4);

[t,q] = ode45(@(t,X) rossler(t,X,a,b,c),tspan,x0,options);

nt = length(t);
tend = floor(nt/2):nt;
figure;
plot(t,q);legend('x','y','z');
xlabel('t');

figure;
plot(q(tend,1),q(tend,2));
xlabel('x');ylabel('y');
