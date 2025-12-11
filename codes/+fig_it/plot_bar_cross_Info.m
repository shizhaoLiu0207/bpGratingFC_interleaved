function [stats_info_cardinal, stats_info_oblique] = plot_bar_cross_Info(results_all, session_list_plot, plotOptions)
global bpGlobal ftsize
idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
    ismember({results_all(:).sessionStr}, session_list_plot);


I_cardinal_real = cell2mat({results_all(idx).fisher_cardinal_cardinal_median});
I_cardinal_cross = cell2mat({results_all(idx).fisher_cardinal_oblique_median});
I_cardinal_shuffle = cell2mat({results_all(idx).fisher_cardinal_cardinal_shuffle_median});
I_cardinal_cross_shuffle = cell2mat({results_all(idx).fisher_cardinal_oblique_shuffle_median});


I_oblique_real = cell2mat({results_all(idx).fisher_oblique_oblique_median});
I_oblique_cross = cell2mat({results_all(idx).fisher_oblique_cardinal_median});
I_oblique_shuffle = cell2mat({results_all(idx).fisher_oblique_oblique_shuffle_median});
I_oblique_cross_shuffle = cell2mat({results_all(idx).fisher_oblique_cardinal_shuffle_median});


I_cardinal_real(isnan(I_cardinal_real)) = [];
I_cardinal_cross(isnan(I_cardinal_cross)) = [];
I_cardinal_shuffle(isnan(I_cardinal_shuffle)) = [];
I_cardinal_cross_shuffle(isnan(I_cardinal_cross_shuffle)) = [];

I_oblique_real(isnan(I_oblique_real)) = [];
I_oblique_cross(isnan(I_oblique_cross)) = [];
I_oblique_shuffle(isnan(I_oblique_shuffle)) = [];
I_oblique_cross_shuffle(isnan(I_oblique_cross_shuffle)) = [];

switch plotOptions.plot_data
    case 'real'
        data_plot_cardinal_real = I_cardinal_real;
        data_plot_cardinal_cross = I_cardinal_cross;
        data_plot_oblique_real  = I_oblique_real;
        data_plot_oblique_cross  = I_oblique_cross;
    case 'shuffle'
        data_plot_cardinal_real = I_cardinal_shuffle;
        data_plot_cardinal_cross = I_cardinal_cross_shuffle;
        data_plot_oblique_real  = I_oblique_shuffle;
        data_plot_oblique_cross = I_oblique_cross_shuffle;
end

hold on;
bar(0.5, mean(data_plot_cardinal_real), 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_cardinal);
bar(2.5, mean(data_plot_cardinal_cross), 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_oblique, 'linewidth',3);

bar(5.5, mean(data_plot_oblique_real), 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_oblique);
bar(7.5, mean(data_plot_oblique_cross), 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_cardinal, 'linewidth',3);

ymean = [mean(data_plot_cardinal_real), mean(data_plot_cardinal_cross), mean(data_plot_oblique_real), mean(data_plot_oblique_cross)];




switch plotOptions.errorbar
    case 'CI_sample'
        cardinal_real_CI            = results_all(idx).fisher_cardinal_cardinal_CI;
        cardinal_cross_CI           = results_all(idx).fisher_cardinal_oblique_CI;
        cardinal_shuffle_CI         = results_all(idx).fisher_cardinal_cardinal_shuffle_CI;
        cardinal_cross_shuffle_CI   = results_all(idx).fisher_oblique_cardinal_shuffle_CI;

        oblique_real_CI            = results_all(idx).fisher_oblique_oblique_CI;
        oblique_cross_CI           = results_all(idx).fisher_oblique_cardinal_CI;
        oblique_shuffle_CI         = results_all(idx).fisher_oblique_oblique_shuffle_CI;
        oblique_cross_shuffle_CI   = results_all(idx).fisher_cardinal_oblique_shuffle_CI;

        switch plotOptions.plot_data
            case 'real'
                y_errorbar_low = ymean - [cardinal_real_CI(1), cardinal_cross_CI(1), ...
                    oblique_real_CI(1), oblique_cross_CI(1)];
                y_errorbar_high  = [cardinal_real_CI(2), cardinal_cross_CI(2),...
                    oblique_real_CI(2), oblique_cross_CI(2)] - ymean;
            case 'shuffle'
                y_errorbar_low = ymean - [cardinal_shuffle_CI(1), cardinal_cross_shuffle_CI(1), ...
                    oblique_shuffle_CI(1), oblique_cross_shuffle_CI(1)];
                y_errorbar_high  = [cardinal_shuffle_CI(2), cardinal_cross_shuffle_CI(2),...
                    oblique_shuffle_CI(2), oblique_cross_shuffle_CI(2)] - ymean;
        end
        errorbar( [0.5,2.5, 5.5,7.5], ymean, y_errorbar_low, y_errorbar_high, '.', 'LineWidth', 2, 'color', 'black');

    case 'SEM_session'

        ysem  = [std(data_plot_cardinal_real), std(data_plot_cardinal_cross),  ...
            std(data_plot_oblique_real), std(data_plot_oblique_cross)] / sqrt(sum(idx));

        errorbar( [0.5,2.5, 5.5,7.5], ymean, ysem, '.', 'LineWidth', 2, 'color', 'black');



end



if plotOptions.dottest
    stats_info_cardinal = fig.show_ttest(data_plot_cardinal_real, data_plot_cardinal_cross, [0.5,2.5]);
    stats_info_oblique  = fig.show_ttest(data_plot_oblique_real, data_plot_oblique_cross, [5.5,7.5]);
    %fig.show_ttest(I_cardinal_real, I_cardinal_cross_corr, [1,3]);


end



set(gca, 'fontsize', plotOptions.ftsize)
set(gca, 'TickLabelInterpreter', 'latex')
switch plotOptions.plot_data
    case 'real'
        set(gca, 'xtick', [0.5,2.5, 5.5,7.5], 'xticklabels', {'$I_\textrm{real}$';'$I_\textrm{cross}$'; ...
            '$I_\textrm{real}$';'$I_\textrm{cross}$'})
    case 'shuffle'
         set(gca, 'xtick', [0.5,2.5, 5.5,7.5], 'xticklabels', {'$I_\textrm{shuffle}$';'$I_\textrm{shuffle,cross}$'; ...
            '$I_\textrm{shuffle}$';'$I_\textrm{shuffle,cross}$'})
end


ylabel('Linear Fisher information','Interpreter','latex')
end