clc; clear; close all;

% Criar grade de pontos
[x, y] = meshgrid(linspace(-2, 2, 50), linspace(-2, 2, 50));

% Criar figura e configuração do colormap
figure;
hold on;
colormap jet;
caxis([-1, 1]); % Ajustar escala de cores

% Criar gráfico inicial
Z = sin(pi*x) .* cos(pi*y);
contourPlot = contourf(x, y, Z, 20);
colorbar;

% Loop para animação
for t = 1:50
    Z = sin(pi*x + 0.1*t) .* cos(pi*y + 0.1*t); % Atualizar os dados
    
    % Remover contornos antigos corretamente
    delete(findall(gca, 'Type', 'Contour'));

    % Criar novos contornos
    contourPlot = contourf(x, y, Z, 20);
    
    pause(0.1); % Pequeno atraso para visualizar a animação
end

hold off;
