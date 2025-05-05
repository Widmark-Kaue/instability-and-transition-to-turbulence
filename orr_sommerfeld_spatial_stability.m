%% Exercise #6: Spatial Stability with tanh profile
clc; clear all; close all

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
plot(real(lambda), imag(lambda), 'ro', 'LineWidth',1.5)
hold off
grid on
xlabel('\alpha_r')
ylabel('\alpha_i')
xlim([-1, 1])
ylim([-2.5, 2.5])
