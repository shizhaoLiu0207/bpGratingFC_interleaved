function stats_info = show_ttest(group1, group2, xplot)


% 
% Perform two-sample t-test
[~, p_value,~,stats] = ttest(group1, group2);
% 
stats_info.mu_1 = mean(group1);
stats_info.mu_2 = mean(group2);
stats_info.std_1 = std(group1);
stats_info.std_2 = std(group2);
stats_info.p_value = p_value;
stats_info.df = stats.df;
stats_info.tstat = stats.tstat;
stats_info.type = 'One-sample t-test';
% Calculate means of each group
mean1 = mean(group1);
mean2 = mean(group2);
% 
sem1 = std(group1) / sqrt(numel(group1));
sem2 = std(group2) / sqrt(numel(group2));
% 
% Plot vertical lines from each bar
y = 1.05 * [mean1 + sem1 , mean2 + sem2];
offset = 0.1 * max(y);  % The vertical line's height above the bar
% 
% Vertical lines from each bar
plot([xplot(1), xplot(1)], [y(1), max(y) + offset], '-k', 'LineWidth', 1.5);  % Vertical line from bar 1
plot([xplot(2), xplot(2)], [y(2), max(y) + offset], '-k', 'LineWidth', 1.5);  % Vertical line from bar 2
% 
% Horizontal line connecting the two vertical lines
plot([xplot(1), xplot(2)], [max(y) + offset, max(y) + offset], '-k', 'LineWidth', 1.5);  % Horizontal line
% 
% % Vertical line connecting the horizontal line to the p-value text
x_mid = mean(xplot);  % Midpoint between the two bars
y_max = max(y) + offset + 0.1 * max(y);  % Position for the p-value text
plot([x_mid, x_mid], [max(y) + offset, y_max], '-k', 'LineWidth', 1.5);  % Vertical line to p-value
% 
% Annotate p-value on top of the vertical line
%text(x_mid, 1.05* y_max, sprintf('p = %.e', p_value), 'HorizontalAlignment', 'center','fontsize', 18,'Interpreter','latex');
% 
text(x_mid, 1.05* y_max,fig_it.p2star(p_value), 'HorizontalAlignment', 'center','fontsize', 18,'Interpreter','latex');
end