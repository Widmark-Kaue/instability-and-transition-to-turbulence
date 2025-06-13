clear all; close all;
%% load base flow
currfile = dir('BF2D_uv_L6p6_long.dat')';
    fname=currfile.name
    data = load(fname);
    x = data(:,1);
    y = data(:,2);
    u = data(:,3);
    v = data(:,4);

    xmin = min(x(:));
    xmax = max(x(:));x
    ymin = min(y(:));
    ymax = 30 ; %   max(y(:));
    
    Nx = 300;
    Ny = 40;

    X = linspace(xmin,xmax,Nx);
%    Y = linspace(ymin,ymax,Ny);
    Y = linspace(0,sqrt(ymax),Ny).^2;

    [yy,xx] = meshgrid(Y ,X);
    U=xx*nan; V=U;
    U(:) = (griddata(x,y,u,xx(:),yy(:)));
    V(:) = (griddata(x,y,v,xx(:),yy(:)));
    
    % enforce periodicity of U and V in x
    Ny = length(Y); Nx = length(X);
    Nper = Nx/10;
    U1 = U(Nx-Nper,:);
    U2 = U(1,:);
    V1 = V(Nx-Nper,:);
    V2 = V(1,:);
    x1 = X(Nx-Nper);
    Lx1x2 = X(end)-x1;

    U(Nx-Nper+1:Nx,:) = (U1.'+(U2-U1).'*(X(Nx-Nper+1:Nx)-x1)/Lx1x2).';
    V(Nx-Nper+1:Nx,:) = (V1.'+(V2-V1).'*(X(Nx-Nper+1:Nx)-x1)/Lx1x2).';


    figure
        subplot(211);contourf(X,Y,U.',128,'linecolor', 'none');colorbar;
        xlabel('x');     ylabel('y'); title('U');

        subplot(212);contourf(X,Y,V.',128,'linecolor', 'none');colorbar;
        xlabel('x');     ylabel('y'); title('V');
        
%% parameters
Re = 450;
nmodes = 100;
BLforcing = false;

beta = 0.166;
%beta = 0;

smax = 0.07; % maximum damping in fringe - set to zero for no fringe
smax = 0.15;

%% grid and derivatives

disp('Forming operators');
tic;

facy = (y(end)-y(1)) / 2;
Lx = (X(end) + X(2));

[Dy,D2y,~] = DiffMat(Ny,Y,5);
dy = Y(2) - Y(1);
Wy = dy*eye(size(Dy));
iWy = 1/dy*eye(size(Dy));


[~,Dx] = fddif(Nx);
[~,D2x] = fd2dif(Nx);

fac = Lx/(2*pi);
Dx = Dx/fac; D2x = D2x/(fac^2);
dx = X(2) - X(1);
Wx = dx*eye(size(Dx));
iWx = 1/dx*eye(size(Dx));


poswall = find(abs(yy)<1e-8 | abs(yy-max(Y))<1e-8);
poswall = poswall(:);


%% full operators
IIy = speye(size(Dy));
IIx = speye(size(Dx));

DY = kron(Dy,IIx);
D2Y = kron(D2y,IIx);

DX = kron(IIy,Dx);
D2X = kron(IIy,D2x);

W = kron(Wy,IIx)*kron(IIy,Wx);
iW = kron(iWy,IIx)*kron(IIy,iWx);

%% base flow:

UU = U(:); % takes U(nx,ny) and transforms it to UU(nx*ny,1)
UUy = DY*UU;
UUx = DX*UU;

VV = V(:);
VVy = DY*VV;
VVx = DX*VV;


Uy = reshape(UUy,Nx,Ny); % takes UUy(nx*ny,1) and transforms it to Uy(nx,ny)
Ux = reshape(UUx,Nx,Ny);


figure;
contourf(xx,yy,Uy,33,'LineStyle','none'); colorbar;xlabel('x');ylabel('y');title('dU/dy(x,y)');


%% fringe zone

xc = 15*Lx/16;
Lfringe = Lx/32;

%xf = xx<3*Lx/4;
xf = xx<748 & xx>50;
%xC = xx>=x0 & xx<=x1;
xC = xf;

sigmaf = smax*exp(-(xx-xc).^2/Lfringe.^2);

figure;
contourf(xx,yy,sigmaf,33); colorbar;xlabel('x');ylabel('y');title('\sigma(x,y)');
%% resolvent operator

II = speye(size(DY)); ZZ = sparse(size(DY,1),size(DY,2));%ZZ = zeros(size(DY));

[np,~] = size(II); % number of points of full grid vector

IIf = spdiags(xf(:),0,np,np); % physical domain
IIC = spdiags(xC(:),0,np,np); % physical domain

if(BLforcing)
    IIbl = diag(UU<0.98); IIf = IIbl*IIf;
end


smat = sigmaf(:); sMat = spdiags(smat,0,np,np);

Lapl = D2X + D2Y -beta^2*II;
Conv = spdiags(UU,0,np,np)*DX + spdiags(VV,0,np,np)*DY;


L = [Conv+spdiags(UUx,0,np,np)-1/Re*Lapl+sMat spdiags(UUy,0,np,np)                     ZZ DX;
      spdiags(VVx,0,np,np)                    Conv+spdiags(VVy,0,np,np)-1/Re*Lapl+sMat ZZ DY;
        ZZ ZZ Conv-1/Re*Lapl+sMat 1i*beta*II;
        DX DY 1i*beta*II ZZ];

%% boundary conditions
[np,~] = size(II); % number of points of full grid vector
L(poswall,:) = 0;
L(poswall,poswall) = eye(length(poswall)); %u=0;
L(np+poswall,:) = 0;
L(np+poswall,np+poswall) = eye(length(poswall)); %v=0;
L(2*np+poswall,:) = 0;
L(2*np+poswall,2*np+poswall) = eye(length(poswall)); %w=0;
%L(3*np+poswall,:) = 0;
%L(3*np+poswall,3*np+poswall) = eye(length(poswall)); %p=0;




   F = [II ZZ ZZ ZZ;
        ZZ II ZZ ZZ;
        ZZ ZZ II ZZ;
        ZZ ZZ ZZ ZZ];

    F(poswall,:) = 0;
    F(np+poswall,:) = 0;
    F(2*np+poswall,:) = 0;
%    B(3*np+poswall,:) = 0;
    

    toc
    disp('Solving eigenvalue problem');
    tic
        [Vi,lambda] = eigs(L,F,nmodes,'smallestabs');
    toc
    
    omega = diag(lambda)/(1i);
    figure;
    plot(real(omega),imag(omega),'ks');
    xlabel('\omega_r'); ylabel('\omega_i');

    [~,ord] = sort(imag(omega),'descend');
    omega = omega(ord);
    Vi = Vi(:,ord);

        
    nmodes = 4;
    for i=1:nmodes
        f = Vi(:,i);

        u = reshape(f(1:np),Nx,Ny);
        v = reshape(f(np+1:2*np),Nx,Ny);
        w = reshape(f(2*np+1:3*np),Nx,Ny);
        
        figure;

        subplot(3,1,1);
        contourf(xx,yy,real(u),33,'LineStyle','none'); colorbar;
        xlabel('x');ylabel('y');title(['u, \omega=' num2str(omega(i))]);

        subplot(3,1,2);
        contourf(xx,yy,real(v),33,'LineStyle','none'); colorbar;
        xlabel('x');ylabel('y');title(['v, \omega=' num2str(omega(i))]);
        subplot(3,1,3);
        contourf(xx,yy,real(w),33,'LineStyle','none'); colorbar;
        xlabel('x');ylabel('y');title(['w, \omega=' num2str(omega(i))]);
    end 