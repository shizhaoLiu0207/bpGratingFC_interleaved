function make_scatter_match(p1, p2, y_plot, rgb_colors, plotOptions)
p1_name     = plotOptions.p1_name;
p2_name     = plotOptions.p2_name;
mapName     = plotOptions.mapName;
max_min     = plotOptions.max_min;
if ~isfield(plotOptions, 'doContour')
    plotOptions.doContour = false;
end
% cardinal_delta_list = unique(cardinal_delta);
% cardinal_prior_list = unique(cardinal_prior);
p1_list = unique(p1);
p2_list = unique(p2);
nx = numel(p1_list); ny = numel(p2_list);
y_plot_all = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        idx = p1 == p1_list(i) & p2 == p2_list(j);
        y_plot_subset = y_plot(idx);
        rgb_colors_subset = rgb_colors(idx,:);
        switch max_min
            case 'min'
                %%%% the minimum value, good for the distance measurement
                [y_plot_all(i,j), idx_best] = min(y_plot_subset);
            case 'max_abs'
                %%%% the value with maximum absolute value, good for the
                %%%% assymetry measurement
                [~, idx_best] = max(abs(y_plot_subset));
                y_plot_all(i,j) = y_plot_subset(idx_best);
        end
        plot(p1_list(i), p2_list(j), '.','Markersize',60,'Color',rgb_colors_subset(idx_best,:)); hold on        
    end
end
if plotOptions.doContour
    empirical_value = plotOptions.empirical_value;
    [p1_grid, p2_grid] = meshgrid(p1_list, p2_list);

    p1_fine = linspace(min(p1_list), max(p1_list), 50);
    p2_fine = linspace(min(p2_list), max(p2_list), 50);
    [p1_fine_grid, p2_fine_grid] = meshgrid(p1_fine, p2_fine);
    
    y_fine = interp2(p1_grid, p2_grid, y_plot_all', ...
        p1_fine_grid, p2_fine_grid);
    
    contour(p1_fine, p2_fine, y_fine, ...
        [empirical_value empirical_value], ...
        'k', 'LineWidth', 2);

    % ;
    % % Add contour where model output matches empirical value
    % [C, h] = contour(p1_list, p2_list, y_plot_all', ...
    %     [empirical_value empirical_value], ...
    %     'k', 'LineWidth', 2);
    
    %clabel(C, h, 'FontSize', 16, 'Color', 'k');
end
clim([min(y_plot_all(:)), max(y_plot_all(:))]);
colormap(mapName)
colorbar('location','northoutside')
title('MSE of Diff. (%)')
set(gca, 'fontsize', 16);

xlabel(p1_name); ylabel(p2_name); 

box off

end