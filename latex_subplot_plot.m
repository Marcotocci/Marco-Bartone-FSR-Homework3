function g = latex_subplot_plot(x, y1, y2, ...
                                 x_label_latex, ...
                                 y1_label_latex, y2_label_latex, ...
                                 y1_legend_latex, y2_legend_latex, ...
                                 title1_latex, title2_latex, ...
                                 pdf_name)
% Plotta due segnali in subplot con LaTeX e salva in PDF

% Impostazioni LaTeX
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')

lw = 2;

% Crea figura
g = figure('Renderer', 'painters', 'Position', [100 100 900 600]);

% Primo subplot
subplot(2,1,1)
plot(x, y1, 'LineWidth', lw, 'Color', [0.2, 0.2, 0.2]);
ylabel(y1_label_latex)
title(title1_latex)
%legend(y1_legend_latex, 'Location', 'north', 'Orientation', 'horizontal', 'AutoUpdate', 'off')
grid on; box on;
set(gca, 'FontSize', 12);

% Secondo subplot
subplot(2,1,2)
plot(x, y2, 'LineWidth', lw, 'Color', [0.2, 0.2, 0.2]);
xlabel(x_label_latex)
ylabel(y2_label_latex)
title(title2_latex)
%legend(y2_legend_latex, 'Location', 'north', 'Orientation', 'horizontal', 'AutoUpdate', 'off')
grid on; box on;
set(gca, 'FontSize', 12);

% Impostazioni generali figura
set(gcf, 'Color', 'w');
set(g, 'MenuBar', 'none');
set(g, 'ToolBar', 'none');

% Esportazione PDF
exportgraphics(g, pdf_name);

end