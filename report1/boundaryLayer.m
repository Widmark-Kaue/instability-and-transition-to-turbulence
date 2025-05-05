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
ReVec = [10:1:100]*100 ;
alpha = 1;

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
%% EigenSpectrum
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
%% Grid convergence Plot
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

%% Wave number variation
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
%% Wave number variation plot 

fig = figure('Position', [100, 100, 800, 300]);
hold on
plot(real(c(1,:)), imag(c(1,:)), 'bs', 'LineWidth',Lwth, 'MarkerFaceColor', 'b', 'DisplayName',b1Label)
plot(real(c(2,:)), imag(c(2,:)), 'ks', 'LineWidth',Lwth, 'DisplayName',b2Label)
annotation('textarrow', [0.17 0.19], [0.6 0.68]);

hold off

xlabel('c_r')
ylabel('c_i')
grid on
box on
legend

