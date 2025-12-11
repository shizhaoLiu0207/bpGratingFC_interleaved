function plot_bar_errorbar(x, mu, CI, plotColor, edge_color)

bar(x, mu, 'facecolor',plotColor, 'EdgeColor',edge_color, 'linewidth',3);
if numel(CI) == 2
    %%%% confidence interval
    errorbar(x, mu, mu - CI(1), CI(2) - mu,'.', 'LineWidth', 2, 'color', 'black');
elseif  numel(CI) == 1
    %%% sem or std
    errorbar(x,mu, CI,'.', 'LineWidth', 2, 'color', 'black');
end