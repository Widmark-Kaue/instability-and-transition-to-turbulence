close all; clear all; clc
path(path, 'src')

%% Constans
Re = 10000; 
alpha = 1;
poiseuilliFlow = @(y) 1-y.^2;
%% grid
N = 360; % degree of highest Chebyshev polynomial
[~, lambda_N360, L_N360, F_N360] = orrSommerfeld(N, alpha, Re, ...
    poiseuilliFlow);

% Plot Eigenspectrum
figure
hold on
plot(real(lambda_N360), imag(lambda_N360), 'ko', 'LineWidth',1.5)
hold off
ylim([-2, 0.1])
xlim([0, 1])

title('Eigenspectrum, \alpha = 1')
xlabel('c_r')
ylabel('c_i')
grid on
%% Grid N = 180
N = 180; % degree of highest Chebyshev polynomial
[~, lambda_N180, L_N180, F_N180] = orrSommerfeld(N, alpha, Re, ...
    poiseuilliFlow);

%% Plot eingspectrum
figure
hold on
plot(real(lambda_N180), imag(lambda_N180), 'bo', 'LineWidth',1.5)
plot(real(lambda_N360), imag(lambda_N360), 'ko', 'LineWidth',1.5)
hold off
ylim([-0.5, 0.1])
xlim([0, 1])

title('Eigenspectrum, \alpha = 1')
xlabel('c_r')
ylabel('c_i')
grid on

%% Sensibility to numerical errors
figure
hold on
plot(real(lambda_N360), imag(lambda_N360), 'ko', 'MarkerSize',3)
for i=1:20
    L1 = L_N180 + 1e-5*rand(size(L_N180));
    F1 = F_N180 + 1e-5*rand(size(F_N180));
    
    % solve eigenvalue problem
    [~,lambda_1]=eig(L1,F1);
    lambda_1 = diag(lambda_1);
    plot(real(lambda_1), imag(lambda_1), 'ko', 'MarkerSize',3)
end
hold off
ylim([-0.5, 0.1])
xlim([0, 1])

title('Eigenspectrum, \alpha = 1')
xlabel('c_r')
ylabel('c_i')

%% Convergence of the most unstable mode
N = 10:5:100;
c = zeros(size(N));

for i = 1:length(N)
    [~, lambda] = orrSommerfeld(N(i), alpha, Re, poiseuilliFlow);
    
    % sort eigenvalues
    [~,pos] = sort(imag(lambda), 'descend');
    lambda = lambda(pos);
    
    % Get the most unstable mode
    c(i) = lambda(1);
end

fig = figure('Position', [100, 100, 800, 400]);

sgtitle('Convergence of most unstable mode, \alpha = 1, Re = 10000')

subplot(1, 2,1)
plot(N, real(c), 'k')
xlabel('N')
ylabel('c_r for the most unstable mode')
ylim([0.235, 0.241])
xticks(0:20:100)

subplot(1, 2,2)
plot(N, imag(c), 'k')
xlabel('N')
ylabel('c_i for the most unstable mode')
xticks(0:20:100)

%% Growth rate of most unstable mode - For loops
N = 100;
alphaVec = 0.1:0.01:1.5;
ReVec = [10:5:100]*100;
% ReVec = [10:15:100]*100;

c = zeros(length(ReVec), length(alphaVec));
for i = 1:length(ReVec)
    Re = ReVec(i);
    for j = 1:length(alphaVec)
        alpha = alphaVec(j);
        [~, lambda] = orrSommerfeld(N, alpha,Re, poiseuilliFlow);
        [~,pos] = sort(imag(lambda), 'descend');
        lambda = lambda(pos);
        
        % Get the most unstable mode
        c(i, j) = lambda(1);
    end
end
%% Growth rate of most unstable mode - Plot
omega_i = repmat(alphaVec, length(ReVec), 1).*imag(c);
colors = {'k', 'b', 'r','k:','b:'};
figure
hold on
for i = 1:length(ReVec(1:5))
    plot(alphaVec, omega_i(i,:), colors{i}, ...
        'DisplayName',['Re = ', num2str(ReVec(i))], ...
        'LineWidth',1.5)
end
hold off
title('Growth rate of most unstable mode: Effect of Reynolds number')
xlabel('\alpha')
ylabel('\omega_i')
grid on
legend()
%% Growth rate of most unstable mode - Countour plot
vec = linspace(-4e-3,4e-3, 10);
[X, Y] = meshgrid(ReVec, alphaVec);
figure
contourf(X, Y, omega_i',vec, 'LineStyle','-')
colormap('jet')
c = colorbar;
caxis([-4e-3 4e-3])
c.Label.String = '\omega_i';

axis tight 
shading interp

xlabel('Re')
ylabel('\alpha')
xticks([2e3:2e3:10e3])

%% Neutral Curve
pos = find(round(omega_i', 3) == 0);
RePlot = X(pos);
alphaPlot = Y(pos);
[~, pos] = sort(RePlot);
alphaPlot = alphaPlot(pos);
RePlot = RePlot(pos);

% [~, pos] = unique(RePlot);
% alphaPlot = alphaPlot(pos);
% RePlot = RePlot(pos);

figure
plot(RePlot, alphaPlot, 'Color',[0 0 0 0.5])

xlabel('Re')
ylabel('\alpha')
xlim([1000, 10000])
ylim([0.2, 1.4])
%% plot eigenmodes

figure;
for i=1:12
subplot(3,4,i);plot(y,V(:,i));title(['Mode ' int2str(i)]);
end