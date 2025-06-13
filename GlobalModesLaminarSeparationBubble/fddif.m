function [x,Dx] = fddif(Nx)
L = 2*pi;
dx = L/(Nx);
x = 0:dx:(L-dx);
x = x.';

um = ones(Nx,1);
Dx = spdiags(um*[1 -8 0 8 -1], -2:2, Nx, Nx)/(12);
Dx(2,Nx) = 1/12;
Dx(1,Nx-1:Nx) = [1 -8]/12;
Dx(Nx-1,1) = -1/12;
Dx(Nx,1:2) = [8 -1]/12;
%Dx(1,1:5) = [-3 4 -1 0 0]/2;
%Dx(2,1:5) = [-0.5 0 0.5 0 0];
%Dx(Nx-1,Nx-4:Nx) = [0 0 -0.5 0 0.5];
%Dx(Nx,Nx-4:Nx) = [0 0 1 -4 3]/2;

Dx = Dx/dx;


