function z = latex_triple_subplot_plot(x, y1, y2, y3, ...
                                       y1_label_latex, y2_label_latex, y3_label_latex, ...
                                       x_label_latex, ...
                                       title_latex)
% Plotta 3 segnali con subplot verticali e testo LaTeX, e salva in PDF

% Impostazioni LaTeX
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')

lw = 2;

% Crea figura
z = figure('Renderer', 'painters', 'Position', [100 100 900 700]);

% Titolo globale (opzionale)
sgtitle(title_latex, 'FontSize', 18);

% Primo subplot
subplot(3,1,1)
plot(x, y1, 'LineWidth', lw, 'Color', [0, 0.4470, 0.7410]);
ylabel(y1_label_latex)
grid on; box on;
set(gca, 'FontSize', 14);

% Secondo subplot
subplot(3,1,2)
plot(x, y2, 'LineWidth', lw, 'Color', [0, 0.4470, 0.7410]);
ylabel(y2_label_latex)
grid on; box on;
set(gca, 'FontSize', 14);

% Terzo subplot
subplot(3,1,3)
plot(x, y3, 'LineWidth', lw, 'Color', [0, 0.4470, 0.7410]);
ylabel(y3_label_latex)
xlabel(x_label_latex)
grid on; box on;
set(gca, 'FontSize', 14);

% Estetica figura
set(gcf, 'Color', 'w');
set(z, 'MenuBar', 'none');
set(z, 'ToolBar', 'none');

end
