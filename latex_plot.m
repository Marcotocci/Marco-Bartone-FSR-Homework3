function h = latex_plot(x, y, x_label_latex, y_label_latex, legend_latex, plot_title_latex, legend_position, pdf_name)
% Funzione per creare plot elegante con testo in LaTeX ed esportazione PDF

% Impostazioni LaTeX
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')

% Spessore linea
lw = 2;

% Crea figura
h = figure('Renderer', 'painters', 'Position', [100 100 900 350]);

% Plot
plot(x, y, 'k-', 'LineWidth', lw, 'Color', [0.2, 0.2, 0.2]);
grid on;
box on;
hold on;

% Etichette e titolo
xlabel(x_label_latex)
ylabel(y_label_latex)
title(plot_title_latex)
%legend(legend_latex, 'Orientation', 'horizontal', 'AutoUpdate', 'off')

% Estetica
set(gca, 'FontSize', 14);
set(gcf, 'Color', 'w');
set(h, 'MenuBar', 'none');
set(h, 'ToolBar', 'none');

% Limiti automatici
%xlim([x(1) x(end)])
%ylim([min(y) max(y)])

% Posizione interna aggiustata
%set(gca, 'InnerPosition', [0.1400 0.32 0.82 0.55])
annotation('rectangle',[0 0 1 1],'Color','w');

% Esportazione
exportgraphics(h, pdf_name);

end
