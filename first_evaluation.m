clc; clear all; close all
path(path, 'src');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1 - Rayleigh equation with Velocity profile of plane jet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
alpha = 1;
theta = [0.05 0.02];
N = 200;
[~,y] = cheb(N);

th1Label = ['\theta = ', num2str(theta(1))];
th2Label = ['\theta = ', num2str(theta(2))];

%% a) Velocity profile
planeJetFlow = @(y, theta) 0.5*( 1 + tanh( (y+0.5)./(2*theta) ) ) - 0.5*( 1 + tanh( (y-0.5)./(2*theta) ) );

planeJetFlow_th1 = @(y) planeJetFlow(y, theta(1));
planeJetFlow_th2 = @(y) planeJetFlow(y, theta(2));

figure
hold on
plot(planeJetFlow_th1(y), y, 'k', ...
    'DisplayName',['\theta = ', num2str(theta(1))])
plot(planeJetFlow_th2(y), y, 'k--', ...
    'DisplayName',['\theta = ', num2str(theta(2))])
hold off

title('Velocity Profiles')
xlabel('U')
ylabel('y')
legend()
grid on
%% b) Growth Rates and Phase Speed of the Unstable modes

[V_th1, lambda_th1, ~, ~, y, D, D2] = rayleigh(N, alpha, planeJetFlow_th1);
[V_th2, lambda_th2] = rayleigh(N, alpha, planeJetFlow_th2);

% Sort for most unstable eigenvalues
[~, pos] = sort(imag(lambda_th1), 'descend');
lambda_th1 = lambda_th1(pos);
V_th1 = V_th1(:, pos);

[~, pos] = sort(imag(lambda_th2), 'descend');
lambda_th2 = lambda_th2(pos);
V_th2 = V_th2(:, pos);

r_th1 = real(lambda_th1);
i_th1 = imag(lambda_th1);

r_th2 = real(lambda_th2);
i_th2 = imag(lambda_th2);

str_f = @(c) sprintf('%.2f', imag(c));

% Plot EigenSpectrum
figure
hold on
plot(r_th1, i_th1, 'bs', 'DisplayName',th1Label)
plot(r_th2, i_th2, 'ks', 'DisplayName',th2Label)


text(r_th1(1)+0.05, i_th1(1), str_f(lambda_th1(1)), 'Color', 'b');
text(r_th1(2)+0.05, i_th1(2), str_f(lambda_th1(2)), 'Color', 'b');

text(r_th2(1)+0.05, i_th2(1), str_f(lambda_th2(1)), 'Color', 'k');
text(r_th2(2)+0.05, i_th2(2), str_f(lambda_th2(2)), 'Color', 'k');


hold off
title(['Eigenspectrum, \alpha = ', num2str(alpha)])
xlabel('c_r')
ylabel('c_i')
legend()
%% d) Eigenfunctions of Unstable modes

figure
hold on
plot(y, abs(V_th1(:, 1)), 'b', 'DisplayName', [th1Label, ' - 1st Mode'])
plot(y, abs(V_th1(:, 2)), 'b--', 'DisplayName', [th1Label, ' - 2nd Mode'])

plot(y, abs(V_th2(:, 1)), 'k', 'DisplayName', [th2Label, ' - 1st Mode'])
plot(y, abs(V_th2(:, 2)), 'k--', 'DisplayName', [th2Label, ' - 2nd Mode'])
hold off
xlabel('y')
grid on
legend()

