clear all; close all;

set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',16);

figure;

as = 0.05:0.0005:0.35;
b = 0.2;
c = 5.7;

x0 = rand(1,3);
figure;
hold on;

for a = as
    

tspan = 0:0.1:1000;
options = odeset('RelTol',1e-4);

[t,q] = ode45(@(t,X) rossler(t,X,a,b,c),tspan,x0,options);


nt = size(q,1);

tend = floor(nt/2):(nt-2);
f = q(tend,:);
compon = 2; vplot = 1;
poscut = 0.0; 
pos1 = find((q(tend,compon)<poscut)&q(tend+1,compon)>poscut);
pos1 = pos1(1:(end-1));
%vec = (f(pos1)+f(pos1+1))/2;
vec = f(pos1,vplot) + (poscut-f(pos1,compon))./(f(pos1+1,compon)-f(pos1,compon)) .* (f(pos1+1,vplot)-f(pos1,vplot));

%plot(Rey*ones(size(pos1)),f(pos1,1),'k.');
plot(a*ones(size(pos1)),vec,'k.');
xlabel('c');ylabel('y');title('Poincare section');
drawnow;
x0 = q(end,:);
a
length(pos1)
end


%plot3(q(floor(nt/2):end,1),q(floor(nt/2):end,2),q(floor(nt/2):end,3));

