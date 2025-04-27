function ffpp=findfpp(x)
global beta
f0 = [0 0 x];
[xout,yout]=ode45(@blasius,[0 10],f0);
ffpp = (yout(size(yout,1),2)-1.)^2;