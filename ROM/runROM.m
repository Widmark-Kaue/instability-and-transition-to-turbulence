close all; clear all;
set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',14);

dovideo = false;

disp('Loading data');
tic;
load('NagataUB_ROM_300.mat');
toc

q0 = qUB; % for upper branch solution at Re=140
%q0 = qLB; % for lower branch solution at Re=140

F = zeros(size(q0));

Re = 140;

L = L0 + Re*L2;


tspan = 0:0.5:5000;

options = odeset('RelTol',1e-6);
disp('Simulation');
tic;
[t,q] = ode45(@(t,X) galerkinsys(t,X,L,QQ,F/Re,Re),tspan,q0,options);
toc

figure;
plot(t,q);
drawnow;

nt = length(t);
tend = floor(nt/2):nt;
figure;
plot(q(tend,6),q(tend,7));
xlabel('$a_6$');
ylabel('$a_7$');

%% video of a slice
if(dovideo)
nmodes = size(L,1);
qfield = zeros(size(phi,1),1,size(phi,3),size(phi,4),nt);
for j=1:nt
qfield(:,:,:,:,j) = phi0(:,1,:,:); % laminar solution
end
for i=1:nmodes    
    disp(['Preparing video, ' int2str(i) ' / 300']);
    for j=1:10:nt
        qfield(:,1,:,:,j) = qfield(:,1,:,:,j)+(phi(:,1,:,:,i))*q(j,i);
    end
end

z = squeeze(Z(1,1,:));
y = squeeze(Y(:,1,1));
x = squeeze(X(1,:,1));

% qfield is ny x nx x nz
figure;
for j=1:10:nt
contourf(z,y,squeeze(qfield(:,1,:,1,j)),64,'LineStyle','none'); axis equal;colormap(bluewhitered);
hold on;
quiver(z,y,squeeze(qfield(:,1,:,3,j)),squeeze(qfield(:,1,:,2,j)),'off');
xlabel('z');ylabel('y');
hold off;
title(['t  = ' num2str(t(j))]);
drawnow;
end

end

