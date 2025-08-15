function stats_info = plot_diff_errorbar(results_all, session_list_plot, plotOptions)

global bpGlobal
idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
    ismember({results_all(:).sessionStr}, session_list_plot);


diff_info_percent_cardinal  = cell2mat({results_all(idx).diff_info_percent_cardinal_median});
diff_info_percent_oblique   = cell2mat({results_all(idx).diff_info_percent_oblique_median});

diff_delta_percent_cardinal =  cell2mat({results_all(idx).diff_delta_percent_cardinal_median});
diff_delta_percent_oblique  =  cell2mat({results_all(idx).diff_delta_percent_oblique_median});


if strcmp(plotOptions.errorbar,'CI_sample')
    diff_info_percent_cardinal_CI  = results_all(idx).diff_info_percent_cardinal_CI;
    diff_info_percent_oblique_CI  = results_all(idx).diff_info_percent_oblique_CI;

    diff_delta_percent_cardinal_CI =  results_all(idx).diff_delta_percent_cardinal_CI;
    diff_delta_percent_oblique_CI  =  results_all(idx).diff_delta_percent_oblique_CI;


end
switch plotOptions.plotdata
    case 'info'
        data_plot_cardinal  = diff_info_percent_cardinal;
        data_plot_oblique   = diff_info_percent_oblique;
        if strcmp(plotOptions.errorbar,'CI_sample')
            CI_plot_cardinal = diff_info_percent_cardinal_CI;
            CI_plot_oblique = diff_info_percent_oblique_CI;
        end
    case 'delta'
        data_plot_cardinal  = diff_delta_percent_cardinal;
        data_plot_oblique   = diff_delta_percent_oblique;
        if strcmp(plotOptions.errorbar,'CI_sample')
            CI_plot_cardinal = diff_delta_percent_cardinal_CI;
            CI_plot_oblique = diff_delta_percent_oblique_CI;
        end
end

hold on
switch plotOptions.errorbar
    case 'CI_sample'
        errorbar(1, data_plot_cardinal, data_plot_cardinal - CI_plot_cardinal(1),CI_plot_cardinal(2) - data_plot_cardinal,...
            'o','LineWidth',3,'Color',bpGlobal.color_list.color_cardinal ,'markersize',plotOptions.markersize);
        errorbar(3, data_plot_oblique,  data_plot_oblique - CI_plot_oblique(1),CI_plot_oblique(2) - data_plot_oblique ,...
            'o','LineWidth',3,'Color',bpGlobal.color_list.color_oblique,'markersize',plotOptions.markersize);




    case 'SEM_session'
        errorbar(1, mean(data_plot_cardinal,'omitnan'), std(data_plot_cardinal,'omitnan')/ sqrt(numel(data_plot_cardinal)),...
            'o','LineWidth',3,'Color',bpGlobal.color_list.color_cardinal,'markersize',plotOptions.markersize);
        errorbar(3, mean(data_plot_oblique,'omitnan'), std(data_plot_oblique,'omitnan')/ sqrt(numel(data_plot_oblique)),...
            'o','LineWidth',3,'Color',bpGlobal.color_list.color_oblique,'markersize',plotOptions.markersize);

end

line([0.5,3.5],[0,0],'linestyle','--','color','black','linewidth',2)



set(gca, 'fontsize', plotOptions.ftsize)
set(gca, 'xtick', [1,3], 'xticklabels', {'\color{red}Cardinal';'\color{blue}Oblique'})
set(gca, 'TickLabelInterpreter','tex')
ylabel('Diff.(\%)','Interpreter','latex')

if strcmp(plotOptions.errorbar, 'SEM_session')
    stats_info.mean_cardinal    = mean(data_plot_cardinal,'omitnan');
    stats_info.mean_oblique     = mean(data_plot_oblique,'omitnan');
    stats_info.std_cardinal     = std(data_plot_cardinal,'omitnan');
    stats_info.std_oblique      = std(data_plot_oblique,'omitnan');
else
    stats_info = [];
end
% switch plotOptions.plotdata
%     case 'info'
%         ylabel('Cross - Real (Cross\%)','Interpreter','latex');
%     case 'delta'
%         ylabel('Within - Cross (Redundancy\%)','Interpreter','latex');
% end


end