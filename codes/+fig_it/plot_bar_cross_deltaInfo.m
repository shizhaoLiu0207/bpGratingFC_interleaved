function [stats_info_cardinal, stats_info_oblique] = plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions)
 global bpGlobal ftsize
    idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
        ismember({results_all(:).sessionStr}, session_list_plot);

    
  
    if plotOptions.plotPercent
        deltaI_cardinal =  cell2mat({results_all(idx).delta_percent_cardinal_cardinal_median});
        deltaI_cardinal_cross =  cell2mat({results_all(idx).delta_percent_oblique_cardinal_median});
        deltaI_oblique = cell2mat({results_all(idx).delta_percent_oblique_oblique_median});
        deltaI_oblique_cross = cell2mat({results_all(idx).delta_percent_cardinal_oblique_median});
    else
        deltaI_cardinal =  cell2mat({results_all(idx).delta_cardinal_cardinal_median});
        deltaI_cardinal_cross =  cell2mat({results_all(idx).delta_oblique_cardinal_median});
        deltaI_oblique = cell2mat({results_all(idx).delta_oblique_oblique_median});
        deltaI_oblique_cross = cell2mat({results_all(idx).delta_cardinal_oblique_median});
    end

    idx_nan_cardinal = isnan(deltaI_cardinal) | isnan(deltaI_cardinal_cross);
    idx_nan_oblique  = isnan(deltaI_oblique) | isnan(deltaI_oblique_cross);
    deltaI_cardinal(idx_nan_cardinal) = [];
    deltaI_cardinal_cross(idx_nan_cardinal) = [];
    deltaI_oblique(idx_nan_oblique) = [];
    deltaI_oblique_cross(idx_nan_oblique) = [];


    hold on
    bar(0.5, mean(deltaI_cardinal, 'omitnan'), 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_cardinal);
    bar(2.5, mean(deltaI_cardinal_cross), 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_oblique, 'linewidth',3);

    bar(5.5, mean(deltaI_oblique), 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_oblique);
    bar(7.5, mean(deltaI_oblique_cross), 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_cardinal, 'linewidth',3);
    

  
    switch plotOptions.errorbar
        case 'CI_sample'
            deltaI_cardinal_CI          = results_all(idx).delta_cardinal_cardinal_CI;
            deltaI_cardinal_cross_CI    = results_all(idx).delta_oblique_cardinal_CI;
            deltaI_oblique_CI           = results_all(idx).delta_oblique_oblique_CI;
            deltaI_oblique_cross_CI     = results_all(idx).delta_cardinal_oblique_CI;

            ymean           = [mean(deltaI_cardinal), mean(deltaI_cardinal_cross), mean(deltaI_oblique), mean(deltaI_oblique_cross)];
            y_errorbar_low  =  ymean - [deltaI_cardinal_CI(1), deltaI_cardinal_cross_CI(1), deltaI_oblique_CI(1), deltaI_oblique_cross_CI(1)];
            y_errorbar_high  = [deltaI_cardinal_CI(2), deltaI_cardinal_cross_CI(2), deltaI_oblique_CI(2), deltaI_oblique_cross_CI(2)] - ymean;
            errorbar([0.5,2.5, 5.5,7.5], ymean, y_errorbar_low, y_errorbar_high, '.', 'LineWidth', 2, 'color', 'black');


        case 'SEM_session'
            ymean = [mean(deltaI_cardinal), mean(deltaI_cardinal_cross), mean(deltaI_oblique), mean(deltaI_oblique_cross)];
            ysem = [std(deltaI_cardinal), std(deltaI_cardinal_cross), std(deltaI_oblique), std(deltaI_oblique_cross)] / sqrt(numel(idx));
        
            errorbar([0.5,2.5, 5.5,7.5], ymean, ysem, '.', 'LineWidth', 2, 'color', 'black');
    end

     if plotOptions.dottest
        stats_info_cardinal = fig.show_ttest(deltaI_cardinal, deltaI_cardinal_cross, [0.5,2.5]);
        stats_info_oblique = fig.show_ttest(deltaI_oblique, deltaI_oblique_cross, [5.5,7.5]);
     end
    set(gca, 'fontsize', ftsize)
    set(gca, 'xtick', [0.5,2.5, 5.5,7.5], 'xticklabels', {'Within';'Cross'; 'Within';'Cross'})
    ylabel('$I_\textrm{redundancy}$','Interpreter','latex')

end