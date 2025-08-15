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

    hold on;
    bar(0.5, mean(I_cardinal_real), 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_cardinal);
    bar(2.5, mean(I_cardinal_cross), 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_oblique, 'linewidth',3);
   % bar(3, mean(I_cardinal_cross_corr), 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_oblique, 'linewidth',3 ,'linestyle', '--');
    bar(5.5, mean(I_oblique_real), 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_oblique);
    bar(7.5, mean(I_oblique_cross), 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_cardinal, 'linewidth',3);

    ymean = [mean(I_cardinal_real), mean(I_cardinal_cross), mean(I_oblique_real), mean(I_oblique_cross)];
    
    if plotOptions.plotShuffle
        bar(3, mean(I_cardinal_shuffle), 'facecolor',bpGlobal.color_list.color_cardinal_light, 'EdgeColor',bpGlobal.color_list.color_cardinal_light);
        bar(4, mean(I_cardinal_cross_shuffle), 'facecolor',bpGlobal.color_list.color_cardinal_light, 'EdgeColor',bpGlobal.color_list.color_oblique, 'linewidth',3)
    
        bar(8, mean(I_oblique_shuffle), 'facecolor',bpGlobal.color_list.color_oblique_light, 'EdgeColor',bpGlobal.color_list.color_oblique_light);
        bar(9, mean(I_oblique_cross_shuffle), 'facecolor',bpGlobal.color_list.color_oblique_light, 'EdgeColor',bpGlobal.color_list.color_cardinal, 'linewidth',3);
        ymean = [mean(I_cardinal_real), mean(I_cardinal_cross),  mean(I_cardinal_shuffle), mean(I_cardinal_cross_shuffle), ...
                mean(I_oblique_real), mean(I_oblique_cross),  mean(I_oblique_shuffle), mean(I_oblique_cross_shuffle)];
    end



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
            if plotOptions.plotShuffle
                y_errorbar_low = ymean - [cardinal_real_CI(1), cardinal_cross_CI(1), cardinal_shuffle_CI(1), cardinal_cross_shuffle_CI(1),...
                                          oblique_real_CI(1), oblique_cross_CI(1), oblique_shuffle_CI(1), oblique_cross_shuffle_CI(1)];
                y_errorbar_high  = [cardinal_real_CI(2), cardinal_cross_CI(2), cardinal_shuffle_CI(2), cardinal_cross_shuffle_CI(2),...
                                          oblique_real_CI(2), oblique_cross_CI(2), oblique_shuffle_CI(2), oblique_cross_shuffle_CI(2)] - ymean;
                errorbar([1:4, 6:9], ymean, y_errorbar_low, y_errorbar_high, '.', 'LineWidth', 2, 'color', 'black');
            else
                y_errorbar_low = ymean - [cardinal_real_CI(1), cardinal_cross_CI(1), ...
                                          oblique_real_CI(1), oblique_cross_CI(1)];
                y_errorbar_high  = [cardinal_real_CI(2), cardinal_cross_CI(2),...
                                          oblique_real_CI(2), oblique_cross_CI(2)] - ymean;
                errorbar( [0.5,2.5, 5.5,7.5], ymean, y_errorbar_low, y_errorbar_high, '.', 'LineWidth', 2, 'color', 'black');
            end
        case 'SEM_session'
            if plotOptions.plotShuffle
                ysem  = [std(I_cardinal_real), std(I_cardinal_cross), std(I_cardinal_shuffle), std(I_cardinal_shuffle_cross),  ...
                            std(I_oblique_real), std(I_oblique_cross), std(I_oblique_shuffle), std(I_oblique_shuffle_cross)] / sqrt(sum(idx));
                errorbar([1:4, 6:9], ymean, ysem, '.', 'LineWidth', 2, 'color', 'black');
            else
                ysem  = [std(I_cardinal_real), std(I_cardinal_cross),  ...
                            std(I_oblique_real), std(I_oblique_cross)] / sqrt(sum(idx));
                errorbar( [0.5,2.5, 5.5,7.5], ymean, ysem, '.', 'LineWidth', 2, 'color', 'black');
            end
            

    end

   

    if plotOptions.dottest
        stats_info_cardinal = fig.show_ttest(I_cardinal_real, I_cardinal_cross, [0.5,2.5]);
        stats_info_oblique  = fig.show_ttest(I_oblique_real, I_oblique_cross, [5.5,7.5]);
        %fig.show_ttest(I_cardinal_real, I_cardinal_cross_corr, [1,3]);

        if plotOptions.plotShuffle
            fig.show_ttest(I_cardinal_shuffle, I_cardinal_cross_shuffle, [3,4]);
            fig.show_ttest(I_oblique_shuffle, I_oblique_cross_shuffle, [8,9]);
        end
    end
    
    

    set(gca, 'fontsize', plotOptions.ftsize)
    set(gca, 'TickLabelInterpreter', 'latex')
    % set(gca, 'xtick', [1:5, 7:11], 'xticklabels', {'Real';'Cross-cov';'Cross-corr';'Shuffle';'Cross shuffle'; ...
    %                                              'Real';'Cross-cov';'Cross-corr';'Shuffle';'Cross shuffle';})
    if plotOptions.plotShuffle
        set(gca, 'xtick', [1:4, 6:9], 'xticklabels', {'$I_{real}$';'$I_{cross}$';'I_{shuffle}';'I_{shuffle,cross}'; ...
                                                     'I_{real}';'I_{cross}';'I_{shuffle}';'I_{shuffle,cross}'})
    else
        set(gca, 'xtick', [0.5,2.5, 5.5,7.5], 'xticklabels', {'$I_\textrm{real}$';'$I_\textrm{cross}$'; ...
                                                     '$I_\textrm{real}$';'$I_\textrm{cross}$'})
    end
   
    ylabel('Linear Fisher information','Interpreter','latex')
end