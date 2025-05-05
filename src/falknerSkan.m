function U = falknerSkan(y,b)
% put beta=0 for Blasius
% put beta>0 for favourable pressure gradient
% put beta<0 for adverse pressure gradient
global beta
beta = b;
N = length(y);
options = optimset('TolX',1e-10);

ypp = fminsearch(@findfpp,0.4,options);
f0 = [0 0 ypp];
[xout,yout]=ode45(@blasius,y(N:-1:1),f0);
U = yout(N:-1:1,2);
end