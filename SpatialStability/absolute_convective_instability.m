%%
clear all
close all
clc

path(path, 'src')

%% Constants
H = 4;
N = 200;

n = 2;
Lambda = -1;

%% Monkewitz profile
[D, y] = cheb(N);

D2 = D*D/(H^2);
y = y*H;

U = velocity_monkewitz(y, n, Lambda);
ddU = D2*U;

figure
hold on
plot(U, y, 'b')
plot(ddU, y, 'r')
hold off
grid on
xlabel('U, ddU')
ylabel('y')

U = U(2:N);
ddU = ddU(2:N);
%% Spatial stability for omega = 1 and Re = 20
omega = 1;
Re = 20;

[~, lambda] = orrSommerfeld2(omega, Re, U, ddU, "H",H);

figure
subplot(1, 2, 1)
plot(real(lambda), imag(lambda), 'ko', 'MarkerSize',4)

xlim([-2 2])
ylim([-2 2])
xlabel('\alpha_r')
ylabel('\alpha_i')
title('Eigenspectrum at \omega = 1')

omega = 1.0152+0.053966*1i;
[~, lambda] = orrSommerfeld2(omega, Re, U, ddU, "H",H);
subplot(1, 2, 2)
plot(real(lambda), imag(lambda), 'ko', 'MarkerSize',4)

xlim([-2 2])
ylim([-2 2])
xlabel('\alpha_r')
ylabel('\alpha_i')
title(['Eigenspectrum at \omega = ', num2str(omega)])
%% Looking for a saddle point (Trying change the real part also)
% sizes = [10:-2:4 4:2:10];
omega_r = sort([1.01:1e-3:1.016, 1.0152]);
omega_i = sort([-0.06:0.01:0.06, 0.053966]);
colors = jet(length(omega_i));
markersizes = linspace(2, 8, length(omega_r));
Re = 20;

% matrix_lambda = zeros(4*length(U), length(omega_i), length(omega_r));
n = 100;
matrix_lambda = zeros(n, length(omega_i), length(omega_r));

[II, JJ] = ndgrid(1:length(omega_i), 1:length(omega_r));

OMEGA_R = repmat(omega_r,length(omega_i), 1);
OMEGA_I = repmat(omega_i',1 ,length(omega_r));
OMEGA = OMEGA_R + 1i*OMEGA_I;


for idx = 1:numel(II)
    omega = OMEGA(idx);
    disp(['omega = ' num2str(omega)]);
    [~, lambda] = orrSommerfeld2(omega, Re, U, ddU, "H",H, "useSparse",true, "mode","smallestabs", ...
        "numberOfEigenvalues",n);
    matrix_lambda(:, II(idx), JJ(idx)) = lambda;
end
%% plot variable in imaginary part
figure
hold on
for idx = 1:numel(II)
    i = II(idx);
    j = JJ(idx);
    
    omega_r = real(OMEGA(idx));
    cond = 'off';
    if mod(idx, length(omega_i)) == 1
        cond = 'on';
    end

    plot(real(matrix_lambda(:, i, j)), imag(matrix_lambda(:, i, j)), 'o', ...
            color = colors(i, :), ...
            MarkerSize=markersizes(j), ...
            DisplayName = ['\omega_r = ' num2str(omega_r)], ...
            HandleVisibility = cond)
end
hold off
colormap(colors)
c = colorbar();
caxis([omega_i(1), omega_i(end)])
c.Label.String = '\omega_i';


% title(['Eigenspectrum at \omega = ', num2str(real(omega))])
xlim([-2 2])
ylim([-2 2])
xlabel('\alpha_r')
ylabel('\alpha_i')
legend()
grid on
grid minor
