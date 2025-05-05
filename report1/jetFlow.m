clc; clear; close all
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

% Parallel pool
% parpool('local', 2);  % ou especificar número de workers: parpool('local', 4);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1 - Rayleigh equation with Velocity profile of plane jet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
alphaVec = 0.1:0.1:1.5;
theta = [0.05 0.02];
N = 200;
H = 1;
[~,y] = cheb(N-1);
y = y*H;

th1Label = ['\theta = ', num2str(theta(1))];
th2Label = ['\theta = ', num2str(theta(2))];

% Velocity Profile
planeJetFlow = @(y, theta) 0.5*( 1 + tanh( (y+0.5)./(2*theta) ) ) - 0.5*( 1 + tanh( (y-0.5)./(2*theta) ) );

planeJetFlow_th1 = @(y) planeJetFlow(y, theta(1));
planeJetFlow_th2 = @(y) planeJetFlow(y, theta(2));

%% a) Velocity profile

fig = figure;
hold on
plot(planeJetFlow_th1(y), y, 'k', 'LineWidth', Lwth,...
    'DisplayName',['\theta = ', num2str(theta(1))])
plot(planeJetFlow_th2(y), y, 'k--', 'LineWidth', Lwth, ...
    'DisplayName',['\theta = ', num2str(theta(2))])
hold off

% title('Velocity Profiles')
xlabel('U')
ylabel('y')
ylim([-1, 1])
legend()
grid on

save_fig(fig, 'jet_velocity_profile.eps')

%% Convergence Grid Test
Nvec = 10:10:N;
c = zeros(4, length(Nvec));
alpha = 1;

for i=1:length(Nvec)
    n = Nvec(i);
    [~, lambda] = rayleigh(n, alpha, planeJetFlow_th1, H);
    [~, pos] = sort(imag(lambda), 'descend');
    lambda = lambda(pos);

    % c(1:2, i) = lambda(1:2);

    [~, lambda2] = rayleigh(n, alpha, planeJetFlow_th2, H);
    [~, pos] = sort(imag(lambda2), 'descend');
    lambda2 = lambda2(pos);
    
    % c(3:4, i) = lambda(1:2);
    c(:, i) = [lambda(1:2)' lambda2(1:2)'];
end
%% Plot Convergende test
fig  = figure('Position', [150, 150, 1000, 300]);
% sgtitle(sprintf('Test of Convergence Grid, \\alpha = %.1f', alpha))
subplot(1,2, 1)
hold on
plot(Nvec, imag(c(1,:)), 'b', 'DisplayName',[th1Label,' - 1st Mode'], 'LineWidth',Lwth)
plot(Nvec, imag(c(2,:)), 'b--', 'DisplayName',[th1Label,' - 2nd Mode'], 'LineWidth',Lwth)

plot(Nvec, imag(c(3,:)), 'k', 'DisplayName',[th2Label,' - 1st Mode'], 'LineWidth',Lwth)
plot(Nvec, imag(c(4,:)), 'k--', 'DisplayName',[th2Label,' - 2nd Mode'], 'LineWidth',Lwth)

hold off
xlabel('N')
ylabel('c_i')
legend('Location','northeast')
grid on

subplot(1,2, 2)
hold on
plot(Nvec, real(c(1,:)), 'b', 'LineWidth',Lwth)
plot(Nvec, real(c(2,:)), 'b--', 'LineWidth',Lwth)

plot(Nvec, real(c(3,:)), 'k',  'LineWidth',Lwth)
plot(Nvec, real(c(4,:)), 'k--', 'LineWidth',Lwth)

hold off
xlabel('N')
ylabel('c_r')
grid on

save_fig(fig, 'grid_convergence_jetFlow_alpha1.eps')

%% b) Growth Rates and Phase Speed of the Unstable modes
% lines 1:2 -> theta 1
% lines 3:4 -> theta 2
c = zeros(4, length(alphaVec));


for i=1:length(alphaVec)
    alpha = alphaVec(i);
    [~, lambda_th1] = rayleigh(N, alpha, planeJetFlow_th1, H);
    [~, lambda_th2] = rayleigh(N, alpha, planeJetFlow_th2, H);
    
    % Sort for most unstable eigenvalues
    [~, pos] = sort(imag(lambda_th1), 'descend');
    lambda_th1 = lambda_th1(pos);
    c(1:2, i) = lambda_th1(1:2);
    if alpha > 1.5
        c(1:2, i) = lambda_th1(2:-1:1);
    end

    [~, pos] = sort(imag(lambda_th2), 'descend');
    lambda_th2 = lambda_th2(pos);
    c(3:4, i) = lambda_th2(1:2);
    if alpha > 2.3
        c(3:4, i) = lambda_th2(2:-1:1);
    end
end

omega_i = repmat(alphaVec, 4, 1).*imag(c);
cr = real(c);
pos = find(omega_i == 0);
if ~isempty(pos)
    omega_i(pos) = NaN;
end

fig = figure('Position', [50, 50, 1000, 600]);

subplot(2, 2, 1)
hold on
plot(alphaVec, omega_i(1, :), 'b', 'LineWidth',Lwth, ...
    'DisplayName', [th1Label, ' - 1st Mode'] )
plot(alphaVec, omega_i(2, :), 'b--','LineWidth',Lwth, ...
    'DisplayName', [th1Label, ' - 2nd Mode'])

plot(alphaVec, omega_i(3, :), 'k', 'LineWidth',Lwth, ...
    'DisplayName', [th2Label, ' - 1st Mode'] )
plot(alphaVec, omega_i(4, :), 'k--','LineWidth',Lwth, ...
    'DisplayName', [th2Label, ' - 2nd Mode'])
hold off

title('a)')
xlabel('\alpha')
ylabel('\omega_i')
xlim([0, max(alphaVec)])
grid on
legend('Location','best')


subplot(2, 2, 3)
hold on
plot(alphaVec, cr(1, :), 'b', 'LineWidth',Lwth, ...
    'DisplayName', [th1Label, ' - 1st Mode'] )
plot(alphaVec, cr(2, :), 'b--','LineWidth',Lwth, ...
    'DisplayName', [th1Label, ' - 2nd Mode'])

plot(alphaVec, cr(3, :), 'k','LineWidth',Lwth, ...
    'DisplayName', [th2Label, ' - 1st Mode'] )
plot(alphaVec, cr(4, :), 'k--','LineWidth',Lwth, ...
    'DisplayName', [th2Label, ' - 2nd Mode'])
hold off

title('b)')
xlabel('\alpha')
ylabel('c_r')
xlim([0, max(alphaVec)])
grid on
% legend('Location','best')

subplot(2,2, [2, 4])

x = 0.4;
y = 0.1;

hold on
plot(real(c(1,:)), imag(c(1, :)), 'bs','MarkerFaceColor','b', ...
    'DisplayName',[th1Label , ' - 1st Mode'])
plot(real(c(2,:)), imag(c(2, :)), 'bo', 'MarkerFaceColor','b', ...
    'DisplayName',[th1Label , ' - 2nd Mode'])
% 
% pos = find(alphaVec == 2);
% plot(real(c(2,pos:end)), imag(c(2, pos:end)), 'ro', 'MarkerFaceColor','r')
% plot(real(c(1,pos:end)), imag(c(1, pos:end)), 'rs', 'MarkerFaceColor','r')


plot(real(c(3,:)), imag(c(3, :)), 'ks', 'MarkerFaceColor','k', ...
    'DisplayName',[th2Label , ' - 1st Mode'])
plot(real(c(4,:)), imag(c(4, :)), 'ko', 'MarkerFaceColor','k', ...
    'DisplayName',[th2Label , ' - 2nd Mode'])

% quiver(x, y, dx, dy, 'k', 'LineWidth', 0.5, 'MaxHeadSize', 1.5)
% text(x + dx + 0.01, y + dy, '\alpha increasing', 'FontSize', 12)
annotation('textarrow', [0.6 0.62], [0.3 0.37]);
annotation('textarrow', [0.87 0.87], [0.7 0.63]);
text(x, y, {'\rightarrow', '\alpha increasing'}, ...
    'HorizontalAlignment','center', ...
    'FontSize',9, ...
    'FontWeight','bold', ...
    'VerticalAlignment','middle',...
    'BackgroundColor', [0.8 0.8 0.8], ...   % Cinza claro
     'EdgeColor','k', ...                   % Borda preta
     'LineWidth', 0.5,...                   % Espessura da bora
     'Margin',3)


hold off
% axis tight
title('c)')
xlabel('c_r')
ylabel('c_i')
legend('Location', 'best')
grid on

save_fig(fig, 'growth_rate_phase_speed_unstable_modes_jetFlow.eps')

%% b) Eigenspectrum for one alpha 
alpha = 1;
[V_th1, lambda_th1, ~, ~, y, D, D2] = rayleigh(200, alpha, planeJetFlow_th1, H);
[V_th2, lambda_th2] = rayleigh(200, alpha, planeJetFlow_th2, H);

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

str_f = @(c) sprintf('(%.2f, %.2f)',real(c), imag(c));

% Plot EigenSpectrum
fig = figure;
hold on
plot(r_th1, i_th1, 'bs', 'DisplayName',th1Label, 'LineWidth',Lwth)
plot(r_th2, i_th2, 'ks', 'DisplayName',th2Label, 'LineWidth',Lwth)


text(r_th1(1)+0.05, i_th1(1), str_f(lambda_th1(1)), 'Color', 'b');
text(r_th1(2)+0.05, i_th1(2), str_f(lambda_th1(2)), 'Color', 'b');

text(r_th2(1)+0.05, i_th2(1), str_f(lambda_th2(1)), 'Color', 'k');
text(r_th2(2)+0.05, i_th2(2), str_f(lambda_th2(2)), 'Color', 'k');


hold off
% title(['Eigenspectrum, \alpha = ', num2str(alpha)])
xlabel('c_r')
ylabel('c_i')
legend('Location','best')

save_fig(fig, 'eigenspectrum_alpha1_jetFlow.eps')

%% d) Eigenfunctions of Unstable modes

figure('Position', [150, 150, 1000, 400])
subplot(1,2,1)
hold on
semilogy(y, abs(V_th1(:, 1)), 'b', 'LineWidth',Lwth,'DisplayName', th1Label)
semilogy(y, abs(V_th2(:, 1)), 'k', 'LineWidth',Lwth, 'DisplayName', th2Label)
hold off
xlabel('y')
ylabel('V(y)')
title('1st Mode')
grid on
legend()


subplot(1, 2,2)
hold on
semilogy(y, abs(V_th1(:, 2)), 'b', 'LineWidth',Lwth,'DisplayName', th1Label)
semilogy(y, abs(V_th2(:, 2)), 'k', 'LineWidth',Lwth,'DisplayName', th2Label)
hold off
xlabel('y')
ylabel('V(y)')
title('2nd Mode')
grid on

