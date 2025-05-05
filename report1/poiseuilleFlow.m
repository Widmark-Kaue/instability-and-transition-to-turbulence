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
% 
% parpool('local', 4);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q2 - Rayleigh and Orr-Sommerfeld equation with Poiseuille Flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
N = 200;
alpha = 1;
Re = 7000;
baseFlow = @ (y) 1 - y.^2;

dd2U = -2;
%% a.2) Velocity Profile
[D, y] = cheb(N-1);
figure
hold on
plot(baseFlow(y), y, 'k', 'DisplayName','U');
hold off
xlabel('Velocity')
ylabel('y')
legend()

%% a) Temporal Stability with inviscid Pertubations
Nvec = 100:50:N;
alphaVec = 0.5:0.5:1.5;
plotpos = reshape(1:9, [3,3]);
fig = figure('Position', [50, 50, 1000, 600]);
for i = 1:length(alphaVec)
    for j= 1:length(Nvec)

    [~, lambda] = rayleigh(Nvec(j), alphaVec(j), baseFlow);
    [~, pos] = sort(imag(lambda), 'descend');
    lambda = lambda(pos);
    
    subplot(length(Nvec), length(alphaVec),plotpos(i, j))
    plot(real(lambda), imag(lambda), 'bs')
    % xlabel('c_r')
    % ylabel('c_i')
    title(['N =' num2str(Nvec(j)) ', \alpha = ' num2str(alphaVec(i))])
    end
end

han = axes(gcf, 'Visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

xlabel(han, 'c_r');
ylabel(han, 'c_i');

save_fig(fig, 'stability_poiseuille_inviscid.eps')
%% b) Grid convergence
Nvec = 10:5:N;
c = zeros(size(Nvec));

for i = 1:length(Nvec)
    [~, lambda] = orrSommerfeld(Nvec(i), alpha, Re, baseFlow);
    
    % sort eigenvalues
    [~,pos] = sort(imag(lambda), 'descend');
    lambda = lambda(pos);
    
    % Get the most unstable mode
    c(i) = lambda(1);
end

fig = figure('Position', [100, 100, 800, 300]);

% sgtitle('Convergence of most unstable mode, \alpha = 1, Re = 10000')

subplot(1, 2,1)
plot(Nvec, real(c), 'k')
xlabel('N')
ylabel('c_r for the most unstable mode')
grid on
% ylim([0.235, 0.241])
% xticks(0:20:100)

subplot(1, 2,2)
plot(Nvec, imag(c), 'k')
xlabel('N')
ylabel('c_i for the most unstable mode')
grid on
% xticks(0:20:100)

save_fig(fig, 'grid_convergence_poiseulle_viscou.eps')
%% b) Temporal Stability with viscous pertubations
[~, lambda] = rayleigh(N, alpha, baseFlow);
[~, lambda2] = orrSommerfeld(N, alpha, Re, baseFlow);
[~, pos] = sort(imag(lambda2), 'descend');
lambda2 = lambda2(pos);


fig = figure;
dx = -0.2;
dy = -0.1;
hold on
plot(real(lambda), imag(lambda), 'ks', ...
    'MarkerSize',3,...
    'DisplayName', 'Inviscid')
plot(real(lambda2), imag(lambda2), 'bs', 'DisplayName','Viscous')
plot(real(lambda2(1)), imag(lambda2(1)), 'bs', 'MarkerFaceColor','b', ...
    'HandleVisibility','off');
hold off

ylim([-2, 0.1])
xlim([0, 1])

xlabel('c_r')
ylabel('c_i')

legend('Location','southeast')
box on
grid on

save_fig(fig, 'comp_eigenspectruns_poiseulle.eps')

%% Values to get the neutral curve
N = 100;
alphaVec = [0.1:0.01:1.5];
ReVec = [10:1:100]*100 ;
% ReVec = [10:15:100]*100;

%% Compute single
c = zeros(length(ReVec), length(alphaVec));

tic
for i = 1:length(ReVec)
    Re = ReVec(i);
    for j = 1:length(alphaVec)
        alpha = alphaVec(j);
        [~, lambda] = orrSommerfeld(N, alpha,Re, baseFlow);
        [~,pos] = sort(imag(lambda), 'descend');
        lambda = lambda(pos);
        
        % Get the most unstable mode
        c(i, j) = lambda(1);
    end
end
elapsedTimeSingle = toc
%% Compute parallel 
c2 = zeros(length(ReVec), length(alphaVec));

[II, JJ] = ndgrid(1:length(ReVec), 1:length(alphaVec));
II = II(:);
JJ = JJ(:);
ReVEC = repmat(ReVec, 1, length(alphaVec));
AlphaVEC = repmat(alphaVec, length(ReVec), 1);

tic
parfor idx = 1:numel(II)
    Re = ReVEC(idx);
    alpha = AlphaVEC(idx);
    [~, lambda] = orrSommerfeld(N, alpha,Re, baseFlow);
    [~,pos] = sort(imag(lambda), 'descend');
    lambda = lambda(pos);
    
    % Get the most unstable mode
    c2(idx) = lambda(1);

end
elapsedTimePar = toc
%% Countor plot
omega_i = repmat(alphaVec, length(ReVec), 1).*imag(c2);
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
%% b.b) Neutral Curve
omega_i = repmat(alphaVec, length(ReVec), 1).*imag(c2);

pos = find(round(omega_i', 4) == 0);
Recr = min(X(pos));

fig = figure;
contourf(X, Y, omega_i',[0 0], 'LineStyle','-','LineWidth',1.5, ...
    'FaceColor','none','color',[.5 .5 .5])

xline(Recr, 'LineWidth',1.5, 'Color','k', 'LineStyle','--')
axis tight 
shading interp

text(Recr-250,0.5, ['Re_{cr} = ', num2str(Recr)], ...
    'Rotation', 90, ...
    'FontWeight','bold', ...
    'FontSize',11)
% title('Neutral Curve (\omega_i = 0)')
xlabel('Re')
ylabel('\alpha')
xlim([1000, 10000])
% xticks(ReVec(1:10:end))

save_fig(fig, 'neutra_curve_poiseuille.eps')