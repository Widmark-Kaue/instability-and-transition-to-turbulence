close all; clear all;
clc

path(path, 'src')

%% grid and derivatives
N = 200; % number of gridpoints

%% Constants
H = 20; % Stretche domain from [-1:1] to [-H:H]
alpha = 0.3;
mixingLayerFlow = @(y) 0.5*(1+tanh(y/2));

%% Solve Rayleigh equation
[V, lambda, ~, ~, y, D, D2] = rayleigh(N, alpha, mixingLayerFlow, H);

[~, pos] = sort(imag(lambda), 'descend');
lambda = lambda(pos);
V = V(:, pos);

%% Plot U velocity field
figure
U = mixingLayerFlow(y);
plot(U, y, 'k');
% ylim([-2,2])
xlabel('U')
ylabel('y')

%% Plot EigenSpectrum
figure
plot(real(lambda), imag(lambda), 'bs')

pos = find(imag(lambda) ~= 0);
r_l = real(lambda);
i_l = imag(lambda);
text(r_l(1)-0.1, i_l(1)-0.05, '(K-H) mode', 'Color','b')
text(r_l(end)-0.2, i_l(end)+0.05, 'Conjugate of (K-H)', 'Color','b')

title('Eigenspectrum of mixing layer, \alpha = 0.1')
xlabel('c_r')
ylabel('c_i')


%% Plot one eigenfunction
j = 1;
figure
semilogy(y, abs(V(:, j)), 'k')
% plot(y, abs(V(:, j)), 'k')


xlabel('y')
title('Eigenfunction number 1')

%% Fields
j = 1;
c = lambda(j);
waveLength = 2*pi/alpha;
dx = waveLength/64;
x = 0:dx:waveLength;
[X, Y] = meshgrid(x, y);
Vy = V(:, j)*1e-2;

vxy = @(t) Vy*exp(1i*alpha*(x - c*t));
dvxydx = @(t) 1i*alpha*vxy(t);

%From Continuity
uxy = @(t) -1/(1i*alpha) * (D*Vy) * exp(1i*alpha*(x - c*t));
duxydy = @(t) -1/(1i*alpha) * (D2*Vy) * exp(1i*alpha*(x - c*t));

vorticity_K_H = @(t) (dvxydx(t) - duxydy(t));
vorticity_BaseFLow = repmat(-D*U, 1, length(x));

vorticity = @(t) real(vorticity_K_H(t)) + vorticity_BaseFLow;
%% Contour only for K-H 
figure
contourf(X, Y, real(vorticity_K_H(0)), 60, 'LineStyle','-')
% % contourf(X, Y, vorticity_BaseFLow, 100, 'LineStyle','none')
%% Countour Velocity fields
t = 0;
figure
contourf(X, Y, real(vxy(t)), 20,  'LineStyle','-');
colormap('jet')
title('v_{xy}')
figure
contourf(X, Y, real(uxy(t)), 20,  'LineStyle','-');
colormap('jet')
title('u_{xy}')



%% Initial stages of roll-up shear layer
figure('Position', [150, 150, 850, 450])
ax = axes;

dt = 1;
for t = 0:dt:40
    % delete(findall(gca, 'Type', 'Contour'));
    vec=-10:0.2:0;
    disp('V_K_H:')
    min(min(vorticity(t)))
    disp('V_BF:')
    min(min(vorticity_BaseFLow))
    contourf(ax, X, Y, vorticity(t),vec,'LineStyle','none');
    colorbar
    % caxis(ax, [min(vec), max(vec)]);
    caxis(ax, [-1, 0]);

    xlabel('x')
    ylabel('y')
    title(['Vorticity, t = ', num2str(t), ' s'])
    pause(dt*1e-100)
end
%% look at eigenfunctions
j = 10 ; %% look at jth eigenfunction
close all
figure;
hold on
for j = 1:10
    plot(y,-real(V(:,j)), 'DisplayName',['Ef nº' int2str(j)]);
end
hold off

xlabel('y'); 

% title(['Eigenfunction number ' int2str(j)]);
legend('NumColumns',2);
%% 
alphaVec = 0:0.1:1;
ci = zeros(size(alphaVec));
cr = zeros(size(alphaVec));
U = 0.5*(1+tanh(y));
for i = 1:length(alphaVec)
    alpha = alphaVec(i);

    alpha2 =  alpha^2*eye(size(D2));
    % set eigenvalue problem
    L = diag(U)*(D2 - alpha2) - diag(D2*U) ;
    F = D2 - alpha2;
    % set boundary conditions v=0 at both ends
    L(1,:) = 0; L(1,1) = 1; F(1,:) = 0;
    L(N,:) = 0; L(N,N) = 1; F(N,:) = 0;

    % solve eigenvalue problem
    [V,lambda]=eig(L, F);
    lambda = diag(lambda);
    
    [~, pos] = min(abs(real(lambda) - 0.5));

    ci(i) = abs(imag(lambda(pos)));
    cr(i) = abs(real(lambda(pos)));
   
end

% Michalke Values
alphaM = [1:-0.1:0, 0.4446];
ciM    = [0.0000, 0.0327, 0.0674, 0.1044, 0.1442, 0.1875, 0.2352, 0.2885, 0.3487, 0.4184, 0.5000, 0.2133];

[~, pos] = sort(alphaM);
alphaM = alphaM(pos);
ciM = ciM(pos);

close all
figure
subplot(1,2,1)
hold on
plot(alphaVec, alphaVec.*ci, 'ko')
plot(alphaM, alphaM.*ciM, 'ks-')
hold off
xticks(alphaVec(1:2:end))
xlabel('\alpha')
ylabel('\alpha c_i')

subplot(1, 2,2)
hold on
plot(alphaVec, ci, 'ko', 'DisplayName', 'Code')
plot(alphaM, ciM, 'ks-', 'DisplayName','Michalke')
hold off
xlabel('\alpha')
ylabel('c_i')
legend()
%% Change the velocity profile

figure
plot(tanh(y), y, 0.5*(1+tanh(y)), y)
legend('Tanh(y)', '1/2 *(1+tanh(y)')
ylim([-5, 5])
xlabel('U')
ylabel('y')




