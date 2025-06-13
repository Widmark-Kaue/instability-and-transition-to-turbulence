function [x,D2x] = fd2dif(Nx)
L = 2*pi;
dx = L/(Nx);
x = 0:dx:(L-dx);
x = x.';

um = ones(Nx,1);
D2x = spdiags(um*[-1 16 -30 16 -1], -2:2, Nx, Nx)/12;
D2x(2,Nx) = -1/12;
D2x(1,Nx-1:Nx) = [-1 16]/12;
D2x(Nx-1,1) = -1/12;
D2x(Nx,1:2) = [16 -1]/12;

% D2x(1,1:5) = [2 -5 4 -1 0];
% D2x(2,1:5) = [1 -2 1 0 0];
% D2x(Nx-1,Nx-4:Nx) = [0 0 1 -2 1];
% D2x(Nx,Nx-4:Nx) = [0 -1 4 -5 2];

D2x = D2x/(dx^2);