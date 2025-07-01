%% Exercise #6: Spatial Stability with tanh profile
clc; clear all; close all
path(path, 'src')
%% Constants
omega = 0.11;
Re = 100;
H = 20;

%% Grid 
N = 180;
[D, D2, D4, y] = dirichletNeumannBCs(N);
y = y*H;
D2 = D2/(H^2);
D4 = D4/(H^4);

%% Velocity profile
[D_b, y_b] = cheb(N);

% Streching 
y_b = y_b*H;
D2_b = D_b*D_b;
D2_b = D2_b/(H^2);

U = 0.5*(1 + tanh(y_b/2));
% sech_y = 1./cosh(y);
% Upprime = -2*sech_y.^2 .* tanh(y);
ddU = D2_b*U;

figure
hold on
plot(U, y_b, 'k', 'DisplayName','U')
plot(ddU, y_b, 'b', 'DisplayName','ddU')
hold off
grid on
% xlabel('U, U''''')
ylabel('y')
legend()

U = U(2:N);
ddU = ddU(2:N);

%% Auxiliary functions
II = eye(size(D));
ZZ = zeros(size(D));

F0 = -1i*omega*D2 - 1/Re * D4;
F1 = 1i*diag(U)*D2 - 1i*diag(ddU)*II;
F2 = 1i*omega*II + 2/Re * D2;
F3 = -1i*diag(U)*II;
F4 = -1/Re * II;

L = [
    ZZ II ZZ ZZ;
    ZZ ZZ II ZZ;
    ZZ ZZ ZZ II;
   -F0 -F1 -F2 -F3;
    ];

F = [
    II ZZ ZZ ZZ;
    ZZ II ZZ ZZ;
    ZZ ZZ II ZZ;
    ZZ ZZ ZZ F4;
    ];

%% Evaluate eigenvalues
[V, lambda] = eig(L, F);
lambda = diag(lambda);

%% Plot
figure
hold on
plot(real(lambda), imag(lambda), 'ko', 'LineWidth',0.5)
hold off
grid on
xlabel('\alpha_r')
ylabel('\alpha_i')
xlim([-1, 1])
ylim([-2.5, 2.5])

%% Test function to solve spatial problem
[~, lambda] = orrSommerfeld2(omega, Re, U, ddU, H = H);
figure
hold on
plot(real(lambda), imag(lambda), 'ko', 'LineWidth',0.5)
hold off
grid on
title('Test Function')
xlabel('\alpha_r')
ylabel('\alpha_i')
xlim([-1, 1])
ylim([-2.5, 2.5])

%% Grid convergence

n = [120, 180, 240];
color = ['b', 'k', 'r'];
marker = ['s', 'o', 'x'];


figure
hold on

for i = 1:length(n)
    N = n(i);

    % Velocity profile
    [D_b, y_b] = cheb(N);
    
    % Streching 
    y_b = y_b*H;
    D2_b = D_b*D_b;
    D2_b = D2_b/(H^2);
    
    U = 0.5*(1 + tanh(y_b/2));
    ddU = D2_b*U;

    U = U(2:N);
    ddU = ddU(2:N);

    [~, lambda] = orrSommerfeld2(omega, Re, U, ddU, H = H);
    plot(real(lambda), imag(lambda), [color(i) marker(i)], ...
        'LineWidth',1.5, ...
        'DisplayName', ['N = ' num2str(N)],'MarkerSize',4)

end
hold off
grid on
legend()
title(['Eigenspectrum, \omega =', num2str(omega)])
xlabel('\alpha_r')
ylabel('\alpha_i')
xlim([-1, 1])
ylim([-2.5, 2.5])

%% Plot using phase velocity

N = 240;

% Velocity profile
[D_b, y_b] = cheb(N);

% Streching 
y_b = y_b*H;
D2_b = D_b*D_b;
D2_b = D2_b/(H^2);

U = 0.5*(1 + tanh(y_b/2));
ddU = D2_b*U;

U = U(2:N);
ddU = ddU(2:N);

[V, lambda] = orrSommerfeld2(omega, Re, U, ddU, H = H);
    
Uc = omega./real(lambda);

figure(Position=[100 100 900 400])
subplot(1,2,1)
plot(real(lambda), imag(lambda), 'rs','LineWidth',1.5,'MarkerSize',4)
title(['Eigenspectrum, \omega =', num2str(omega)])
xlabel('\alpha_r')
ylabel('\alpha_i')
xlim([-1, 1])
ylim([-2.5, 2.5])
grid on

subplot(1,2,2)
plot(Uc, imag(lambda), 'rs','LineWidth',1.5,'MarkerSize',4)
title(['Eigenspectrum, \omega =', num2str(omega)])
xlabel('U_c')
ylabel('\alpha_i')
xlim([0, 1])
ylim([-0.2, 1])
grid on
%% Plot Eigenfunction
% Get position of Kelvin-Helmholtz mode
[~, loc] = ismembertol(0.225408, real(lambda), 1e-5);
Vkh = abs(V(1:N-1, loc));
% Vmode = abs(V(1:N-1, 2:4));


% figure(Position=[100 100 1000 400])
figure
% subplot(1, 2, 1)
hold on
plot(U, y_b(2:N), 'k', 'DisplayName', 'Velocity Profile', 'LineWidth',1.5)
plot(Vkh, y_b(2:N), 'b:', 'DisplayName','Kelvin-Helmholtz mode', 'LineWidth',1.5)
hold off
grid on
legend()
xlim([0 1.2])
ylabel('y')
xlabel('U, abs(v)')

% subplot(1, 2, 2)
% hold on
% plot(U, y_b(2:N), 'k', 'DisplayName', 'Velocity Profile', 'LineWidth',1.5)
% plot(Vmode(:, 1), y_b(2:N), 'b:', 'DisplayName','mode 2', 'LineWidth',1.5)
% plot(Vmode(:, 2), y_b(2:N), 'g:', 'DisplayName','mode 3', 'LineWidth',1.5)
% plot(Vmode(:, 3), y_b(2:N), 'r:', 'DisplayName','mode 4', 'LineWidth',1.5)
% 
% hold off
% grid on
% legend()
% xlim([0 1.2])
% ylabel('y')
% xlabel('U, abs(v)')
%% Brigg's criteria

omegai = 0:0.02:0.14;
N = 240;
figure

for i = 1:length(omegai)
    % Velocity profile
    [D_b, y_b] = cheb(N);
    
    % Streching 
    y_b = y_b*H;
    D2_b = D_b*D_b;
    D2_b = D2_b/(H^2);
    
    U = 0.5*(1 + tanh(y_b/2));
    ddU = D2_b*U;
                         
    U = U(2:N);
    ddU = ddU(2:N);
    
    OMEGA = omega + 1i*omegai(i);
    [~, lambda] = orrSommerfeld2(OMEGA, Re, U, ddU, H = H);
    subplot(2, 4, i)
    hold on
    plot(real(lambda), imag(lambda), 'ks', 'LineWidth',1.5,'MarkerSize',4)
    plot(real(lambda(loc)), imag(lambda(loc)), 'ro', 'LineWidth', 1.5, 'MarkerSize', 4)
    hold off
    grid on
    title(['Eigenspectrum, \omega =', num2str(OMEGA)])
    xlabel('\alpha_r')
    ylabel('\alpha_i')
    xlim([-1, 1])
    ylim([-2.5, 2.5])
end
