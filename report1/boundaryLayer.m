clc; 
clear; 
close all
path(path, 'src');

% Path to save images
PATH_SAVE = 'images_first_report';
if ~exist(PATH_SAVE,'dir')
    mkdir(PATH_SAVE)
end
save_fig = @ (fig, name) exportgraphics(fig, ...
        fullfile(PATH_SAVE, name), ...
        'Resolution',500);

% LineWidth
Lwth = 1;

% parpool('local', 4);

% global beta
%% Constants 
alphaVec = [0.1:0.01:1.5];
ReVec = [10:1:50]*100 ;
alpha = 1;
Re = 4500;

N = 300;
H = 20;
h = 1;
[D_bl,y] = cheb(N-1);
y = (y+h)*H;
D_bl = D_bl/H;

b1Label = '\beta = 0.1';
b2Label = '\beta = -0.1';

falknerSkan_b1 = @(y) falknerSkan(y, 0.1);
falknerSkan_b2 = @(y) falknerSkan(y, -0.1);

%% Velocity Profile - beta = 0.1, -0.1
U1 = falknerSkan_b1(y);
U2 = falknerSkan_b2(y);

delta1 = interp1(U1, y, 0.99);
delta2 = interp1(U2, y, 0.99);

fact1 = H/delta1;
fact2 = H/delta2;

% Normalize operators
y_d1 = y/delta1;
y_d2 = y/delta2;
D_d1     = D_bl*delta1;
D_d2     = D_bl*delta2;

% Derivatives
ddU_1 = D_d1*D_d1*U1;
ddU_2 = D_d2*D_d2*U2;


fig = figure('Position', [100, 100, 800, 300]);
subplot(1, 2, 1)
hold on
plot(U1,y_d1,'b', 'LineWidth',Lwth, 'DisplayName','U')
plot(ddU_1, y_d1, 'r', 'LineWidth',Lwth, 'DisplayName', 'U''''')
hold off

xlim([-3, 1.5]);
ylim([0, 10]);

title(b1Label)
xlabel('Velocity')
ylabel('y/\delta')
grid on
box on
legend('Location','west')

subplot(1, 2, 2)
hold on
plot(U2,y_d2,'b', 'LineWidth',Lwth, 'DisplayName','U')
plot(ddU_2, y_d2, 'r', 'LineWidth',Lwth, 'DisplayName', 'U''''')
hold off

xlim([-3, 1.5]);
ylim([0, 10]);


title(b2Label)
xlabel('Velocity')
ylabel('y/\delta')
grid on
box on
% legend('')

save_fig(fig, 'Falkner_Skan_velocity_profile.eps');
%% a) EigenSpectrum
[~, lambda1] =  rayleigh(N, alpha, U1, fact1, h);
[~, pos] = sort(imag(lambda1), 'descend');
lambda1 = lambda1(pos);

[~, lambda2] =  rayleigh(N, alpha, U2, fact2, h);
[~, pos] = sort(imag(lambda2), 'descend');
lambda2 = lambda2(pos);

fig = figure('Position', [100, 100, 800, 300]);
subplot(1, 2, 1)
hold on
plot(real(lambda1), imag(lambda1),'ks', 'DisplayName',b1Label)
hold off
box on

% xlim([0, 1.2])
% ylim([-0.02, 0.02])
xlabel('c_r')
ylabel('c_i')
legend('Location','northwest')

subplot(1, 2, 2)
hold on
plot(real(lambda2), imag(lambda2),'ks', 'DisplayName',b2Label)
hold off
box on
% xlim([0, 1.2])

xlabel('c_r')
ylabel('c_i')
legend('Location','northwest')
save_fig(fig, 'Falkner_Skan_eigenspectrum_alpha1.eps')

%% a) Grid convergence of most unstable mode
Nvec = 10:10:N;
c = zeros(2, length(Nvec));
alpha = 1;

for i=1:length(Nvec)
    [~,y] = cheb(Nvec(i)-1);
    y = (y+h)*H;
    
    U1 = falknerSkan_b1(y);
    U2 = falknerSkan_b2(y);

    delta1 = interp1(U1, y, 0.99);
    delta2 = interp1(U2, y, 0.99);
    
    fact1 = H/delta1;
    fact2 = H/delta2;

    
    [~, lambda1] =  rayleigh(Nvec(i), alpha, U1, fact1, h);
    [~, pos] = sort(imag(lambda1), 'descend');
    lambda1 = lambda1(pos);
    
    [~, lambda2] =  rayleigh(Nvec(i), alpha, U2, fact2, h);
    [~, pos] = sort(imag(lambda2), 'descend');
    lambda2 = lambda2(pos);
    
    c(1:2, i) = [lambda1(1) lambda2(1)];
end 
%% a) Grid convergence Plot
fig = figure('Position', [100, 100, 800, 300]);
subplot(1, 2, 1)
hold on
plot(Nvec, imag(c(1,:)),'k', 'LineWidth',Lwth,'DisplayName',b1Label)
plot(Nvec, imag(c(2,:)),'k--', 'LineWidth',Lwth, 'DisplayName',b2Label)
hold off
xlabel('N')
ylabel('c_i')
legend()
box on
grid on

subplot(1, 2, 2)
hold on
plot(Nvec, real(c(1,:)),'k', 'LineWidth',Lwth,'DisplayName',b1Label)
plot(Nvec, real(c(2,:)),'k--','LineWidth',Lwth, 'DisplayName',b2Label)
hold off
xlabel('N')
ylabel('c_r')
grid on
box on

save_fig(fig, 'Falkner_Skan_grid_Convergence_alpha1.eps')

%% sa) Wave number variation
c = zeros(2, length(alphaVec));

for i=1:length(alphaVec)
    alpha = alphaVec(i);
    [~, lambda1] =  rayleigh(N, alpha, U1, fact1, h);
    [~, pos] = sort(imag(lambda1), 'descend');
    lambda1 = lambda1(pos);
    
    [~, lambda2] =  rayleigh(N, alpha, U2, fact2, h);
    [~, pos] = sort(imag(lambda2), 'descend');
    lambda2 = lambda2(pos);
    
    c(1:2, i) = [lambda1(1) lambda2(1)];
end
%% a) Wave number variation plot 

fig = figure('Position', [100, 100, 800, 300]);
hold on
plot(real(c(1,:)), imag(c(1,:)), 'bs', 'LineWidth',Lwth, 'MarkerFaceColor', 'b', 'DisplayName',b1Label)
plot(real(c(2,:)), imag(c(2,:)), 'ks', 'LineWidth',Lwth, 'DisplayName',b2Label)
annotation('textarrow', [0.18 0.2], [0.6 0.68], 'String',{'\alpha', 'increasing'}, ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold');

annotation('textarrow', [0.27 0.29], [0.68 0.6]);
hold off

xlabel('c_r')
ylabel('c_i')
grid on
box on
legend
save_fig(fig, 'Falkner_Skan_eigenspectrum_variable_alpha.eps')

%% Velocity for Viscou effects
N = 300;
alphaVec = [0.1:0.1:1.5];
ReVec = [10:1:50]*100 ;
[~,y] = cheb(N);
y = (y+h)*H;

U1_visc = falknerSkan_b1(y);
U2_visc = falknerSkan_b2(y);
%% b) Eigenspectrum with viscous effects
[~, lambda1] =  rayleigh(N, alpha, U1, fact1, h);

[~, lambda2] =  rayleigh(N, alpha, U2, fact2, h);

[~, lambda1_visc] =  orrSommerfeld(N, alpha, Re, U1_visc(2:N), fact1, h);

[~, lambda2_visc] =  orrSommerfeld(N, alpha, Re, U2_visc(2:N), fact2, h);


[~, pos] =  sort(imag(lambda1_visc), 'descend');
lambda1_visc = lambda1_visc(pos);
[~, pos] =  sort(imag(lambda2_visc), 'descend');
lambda2_visc = lambda2_visc(pos);
%% b) Eigenspectrum with viscous effects plot

fig = figure('Position', [100, 100, 800, 300]);
subplot(1, 2, 1)
hold on
plot(real(lambda1), imag(lambda1),'ks', ...
    'DisplayName','Rayleigh')
plot(real(lambda1_visc), imag(lambda1_visc), 'bs', ...
    'DisplayName','Orr-Sommerfeld' )
hold off

% xlim([0, 1.2])
ylim([-2, 1])
title(b1Label)
xlabel('c_r')
ylabel('c_i')
box on

legend('Location','northwest')

subplot(1, 2, 2)
hold on
plot(real(lambda2), imag(lambda2),'ks', ...
    'DisplayName','Rayleigh')
plot(real(lambda2_visc), imag(lambda2_visc), 'bs', ...
    'DisplayName','Orr-Sommerfeld' )

plot(real(lambda2_visc(1)), imag(lambda2_visc(1)), 'bs', 'MarkerFaceColor','r')
str = sprintf('(%.2f, %.2f)', real(lambda2_visc(1)), imag(lambda2_visc(1)));
text(real(lambda2_visc(1))-0.15, imag(lambda2_visc(1))+0.2, str, 'FontWeight','bold', 'Color','r')
hold off
% xlim([0, 1.2])
ylim([-2, 1])

title(b2Label)
xlabel('c_r')
ylabel('c_i')
% legend('Location','northwest')
box on

save_fig(fig, 'Falkner_Skan_eigenspectrum_alpha1_Re4500.eps')

%% b) Grid convergence for viscou Falkner Skan

Nvec = 10:10:N;
c = zeros(2, length(Nvec));

for i=1:length(Nvec)
    [~,y] = cheb(Nvec(i));
    y = (y+h)*H;
    
    U1_visc = falknerSkan_b1(y);
    U2_visc = falknerSkan_b2(y);

    delta1 = interp1(U1_visc, y, 0.99);
    delta2 = interp1(U2_visc, y, 0.99);
    
    fact1 = H/delta1;
    fact2 = H/delta2;

    
    [~, lambda1_visc] =  orrSommerfeld(Nvec(i), alpha, Re, U1_visc(2:Nvec(i)), fact1, h);
    [~, pos] =  sort(imag(lambda1_visc), 'descend');
    lambda1_visc = lambda1_visc(pos);

    [~, lambda2_visc] =  orrSommerfeld(Nvec(i), alpha, Re, U2_visc(2:Nvec(i)), fact2, h);
    [~, pos] =  sort(imag(lambda2_visc), 'descend');
    lambda2_visc = lambda2_visc(pos);
    
    c(1:2, i) = [lambda1_visc(1) lambda2_visc(1)];
end 

%% b) Grid convergence Plot - Viscous
fig = figure('Position', [100, 100, 800, 300]);
subplot(1, 2, 1)
hold on
plot(Nvec, imag(c(1,:)),'k', 'LineWidth',Lwth,'DisplayName',b1Label)
plot(Nvec, imag(c(2,:)),'k--', 'LineWidth',Lwth, 'DisplayName',b2Label)
hold off
xlabel('N')
ylabel('c_i')
legend()
box on
grid on

subplot(1, 2, 2)
hold on
plot(Nvec, real(c(1,:)),'k', 'LineWidth',Lwth,'DisplayName',b1Label)
plot(Nvec, real(c(2,:)),'k--','LineWidth',Lwth, 'DisplayName',b2Label)
hold off
xlabel('N')
ylabel('c_r')
grid on
box on

save_fig(fig, 'Falkner_Skan_grid_Convergence_alpha1_Re4500.eps')

%% Compute parallel 
% parpool('local', 4);
disp(' -> Starting Parallel: ')
N = 150;
alphaVec = [0.1:0.01:1.5];
ReVec = [3:1:50]*100 ;
[~,y] = cheb(N);
y = (y+h)*H;

U1_visc = falknerSkan_b1(y);
U2_visc = falknerSkan_b2(y);
c_b1 = zeros(length(ReVec), length(alphaVec));
c_b2 = zeros(length(ReVec), length(alphaVec));

[II, JJ] = ndgrid(1:length(ReVec), 1:length(alphaVec));
II = II(:);
JJ = JJ(:);
ReVEC = repmat(ReVec, 1, length(alphaVec));
AlphaVEC = repmat(alphaVec, length(ReVec), 1);

tic
parfor idx = 1:numel(II)
    Re = ReVEC(idx);
    alpha = AlphaVEC(idx);
    

    [~, lambda1_visc] =  orrSommerfeld(N, alpha, Re, U1_visc(2:N), fact1, h);
    [~, pos] =  sort(imag(lambda1_visc), 'descend');
    lambda1_visc = lambda1_visc(pos);

    [~, lambda2_visc] =  orrSommerfeld(N, alpha, Re, U2_visc(2:N), fact2, h);
    [~, pos] =  sort(imag(lambda2_visc), 'descend');
    lambda2_visc = lambda2_visc(pos);
    
    % [~, lambda] = orrSommerfeld(N, alpha,Re, baseFlow);
    % [~,pos] = sort(imag(lambda), 'descend');
    % lambda = lambda(pos);
    
    % Get the most unstable mode
    c_b1(idx) = lambda1_visc(1);
    c_b2(idx) = lambda2_visc(1);

end
elapsedTimePar = toc

%% Compute single 
% parpool('local', 4);
disp(' -> Starting Single: ')
N = 150;
alphaVec = [0.1:0.05:1.5];
ReVec = [10:1:50]*100 ;
[~,y] = cheb(N);
y = (y+h)*H;

U1_visc = falknerSkan_b1(y);
U2_visc = falknerSkan_b2(y);

c_b1 = zeros(length(ReVec), length(alphaVec));
c_b2 = zeros(length(ReVec), length(alphaVec));

[II, JJ] = ndgrid(1:length(ReVec), 1:length(alphaVec));
II = II(:);
JJ = JJ(:);
ReVEC = repmat(ReVec, 1, length(alphaVec));
AlphaVEC = repmat(alphaVec, length(ReVec), 1);

tic
for idx = 1:numel(II)
    Re = ReVEC(idx);
    alpha = AlphaVEC(idx);
    

    [~, lambda1_visc] =  orrSommerfeld(N, alpha, Re, U1_visc(2:N), fact1, h);
    [~, pos] =  sort(imag(lambda1_visc), 'descend');
    lambda1_visc = lambda1_visc(pos);

    [~, lambda2_visc] =  orrSommerfeld(N, alpha, Re, U2_visc(2:N), fact2, h);
    [~, pos] =  sort(imag(lambda2_visc), 'descend');
    lambda2_visc = lambda2_visc(pos);
    
    % [~, lambda] = orrSommerfeld(N, alpha,Re, baseFlow);
    % [~,pos] = sort(imag(lambda), 'descend');
    % lambda = lambda(pos);
    
    % Get the most unstable mode
    c_b1(idx) = lambda1_visc(1);
    c_b2(idx) = lambda2_visc(1);

end
elapsedTimePar = toc

%% Omega_i
omega_i_b1 = repmat(alphaVec, length(ReVec), 1).*imag(c_b1);
omega_i_b2 = repmat(alphaVec, length(ReVec), 1).*imag(c_b2);

[X, Y] = meshgrid(ReVec, alphaVec);

%% b) Neutral Curve

% pos = find(round(omega_i', 4) == 0);
% Recr = min(X(pos));

fig = figure('Position', [100, 100, 800, 300]);
subplot(1, 2, 1)
hold on
Cneutral = contourf(X, Y, omega_i_b1',[0 0], 'LineStyle','-','LineWidth',1.5, ...
    'FaceColor','none','color',[.5 .5 .5]);
Recr = min(Cneutral(1,2:end));
xline(Recr, 'LineWidth',1.5, 'Color','k', 'LineStyle','--')
axis tight 
shading interp

text(Recr-250,0.5, ['Re_{cr} = ', num2str(int32(Recr))], ...
    'Rotation', 90, ...
    'FontWeight','bold', ...
    'FontSize',9)
title([b1Label])
xlabel('Re')
ylabel('\alpha')
xlim([50, 5000])
box on
% xticks(ReVec(1:10:end))

subplot(1, 2,2)
hold on
Cneutral = contourf(X, Y, omega_i_b2',[0 0], 'LineStyle','-','LineWidth',1.5, ...
    'FaceColor','none','color',[.5 .5 .5]);

Recr = min(Cneutral(1,2:end));
xline(Recr, 'LineWidth',1.5, 'Color','k', 'LineStyle','--')
hold off
axis tight 
shading interp

text(Recr-250,0.5, ['Re_{cr} = ', num2str(int32(Recr))], ...
    'Rotation', 90, ...
    'FontWeight','bold', ...
    'FontSize',9)
title([ b2Label])
xlabel('Re')
ylabel('\alpha')
xlim([50, 5000])
box on
% xticks(ReVec(1:10:end))

save_fig(fig, 'neutra_curve_Falkner_Skan.eps')
