clc 
close all
N = 120;

%% construct Chebyshev differentiation matrix
[D,x]=cheb(N);

%% test derivative
f = exp(x.^2);
df = D*f;
dfan = 2*x.*(exp(x.^2));
plot(x,df,'bs',x,dfan,'k-');
xlabel('x');ylabel('f''');
legend('Numerical derivative','Analytical derivative');

max(abs(df-dfan)) %maximum error

