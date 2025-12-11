function [stats_info,h, data_plot] = plot_diff_errorbar_avg(results_all, session_list_plot, plotOptions)

global bpGlobal
idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
    ismember({results_all(:).sessionStr}, session_list_plot);


diff_info_percent  = cell2mat({results_all(idx).diff_info_percent_median});

diff_delta_percent =  cell2mat({results_all(idx).diff_delta_percent_median});


if strcmp(plotOptions.errorbar,'CI_sample')
    diff_info_percent_CI  = results_all(idx).diff_info_percent_CI;
    diff_delta_percent_CI =  results_all(idx).diff_delta_percent_CI;    
end

switch plotOptions.plotdata
    case 'info'
        data_plot  = diff_info_percent;
        if strcmp(plotOptions.errorbar,'CI_sample')
            CI_plot = diff_info_percent_CI;
           
        end
        plot_color = [131, 201, 117] / 255;
    case 'delta'
        data_plot = diff_delta_percent;
        if strcmp(plotOptions.errorbar,'CI_sample')
            CI_plot = diff_delta_percent_CI;
        end
        plot_color = [208, 201, 124] / 255;
end
switch plotOptions.subjectCode
    case 'Ro'
        plotMarker = bpGlobal.marker.rolo;
    case 'Gr'
        plotMarker = bpGlobal.marker.gremlin;
end
hold on
switch plotOptions.errorbar
    case 'CI_sample'
        errorbar(plotOptions.x_val, data_plot, data_plot - CI_plot(1),CI_plot(2) - data_plot,...
            plotMarker,'LineWidth',3,'Color',plot_color ,'markersize',plotOptions.markersize);

    case 'SEM_session'
        [~,p_val,~,tstats]    = ttest(data_plot);
        h = errorbar(plotOptions.x_val, mean(data_plot,'omitnan'), std(data_plot,'omitnan')/ sqrt(numel(data_plot)),...
            plotMarker,'LineWidth',3,'Color',plot_color,'markersize',plotOptions.markersize);
        % text(plotOptions.x_val, mean(data_plot,'omitnan') + std(data_plot,'omitnan')/ sqrt(numel(data_plot)), ...
        %     fig.p2star(p_val))
end


set(gca, 'fontsize', plotOptions.ftsize)
set(gca, 'TickLabelInterpreter','tex')
ylabel('Diff.(\%)','Interpreter','latex')

if strcmp(plotOptions.errorbar, 'SEM_session')
    
    stats_info.mean_avg   = mean(data_plot,'omitnan');
    stats_info.std_avg    = std(data_plot,'omitnan');
    stats_info.tstats     = tstats.tstat;
    stats_info.df         = tstats.df;
    stats_info.p_val      = p_val;
else
    stats_info = [];
end


end