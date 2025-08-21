function  plot_diff_withintrial(results_all, session_list_plot, plotOptions)

global bpGlobal
nTimebin = 8;
[data_plot_cardinal,data_plot_oblique,CI_plot_cardinal,CI_plot_oblique ] = deal(cell(nTimebin,1));
for t = 1:nTimebin
    idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == t &...
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
            data_plot_cardinal{t}  = diff_info_percent_cardinal;
            data_plot_oblique{t}   = diff_info_percent_oblique;
            if strcmp(plotOptions.errorbar,'CI_sample')
                CI_plot_cardinal{t} = diff_info_percent_cardinal_CI;
                CI_plot_oblique{t} = diff_info_percent_oblique_CI;
            end
            ylabel_str = 'Diff. Information(\%)';
        case 'delta'
            data_plot_cardinal{t}  = diff_delta_percent_cardinal;
            data_plot_oblique{t}   = diff_delta_percent_oblique;
            if strcmp(plotOptions.errorbar,'CI_sample')
                CI_plot_cardinal{t} = diff_delta_percent_cardinal_CI;
                CI_plot_oblique{t} = diff_delta_percent_oblique_CI;
            end
            ylabel_str = 'Diff. Redundancy(\%)';
    end


end
nSession_cardinal       = cellfun(@numel, data_plot_cardinal);
nSession_oblique        = cellfun(@numel, data_plot_oblique);
mean_cardinal           = cellfun(@mean, data_plot_cardinal);
mean_oblique            = cellfun(@mean, data_plot_oblique);
SEM_cardinal            = cellfun(@std, data_plot_cardinal) ./ sqrt(nSession_cardinal);
SEM_oblique             = cellfun(@std, data_plot_oblique) ./ sqrt(nSession_oblique);
if strcmp(plotOptions.errorbar,'CI_sample')
    CI_cardinal             = cat(1,CI_plot_cardinal{:});
    CI_oblique              = cat(1,CI_plot_oblique{:});
end


line([0.5,nTimebin + 0.5], [0,0],'color','black','linestyle','--'); hold on
switch plotOptions.task
    case 'cardinal'
        switch plotOptions.errorbar
            case 'CI_sample'
               h = errorbar([1:nTimebin], mean_cardinal, mean_cardinal - CI_cardinal(:,1),  CI_cardinal(:,2) - mean_cardinal,...
                    'LineWidth',3,'Color',bpGlobal.color_list.color_cardinal);
            case 'SEM_session'
                h = errorbar([1:nTimebin], mean_cardinal, SEM_cardinal,...
                    'LineWidth',3,'Color',bpGlobal.color_list.color_cardinal);
        end
        %egend(h,'Cardinal task')
    case 'oblique'
        switch plotOptions.errorbar
            case 'CI_sample'
                h = errorbar([1:nTimebin], mean_oblique, mean_oblique - CI_oblique(:,1),  CI_oblique(:,2) - mean_oblique,...
                    'LineWidth',3,'Color',bpGlobal.color_list.color_oblique);
            case 'SEM_session'
                h = errorbar([1:nTimebin], mean_oblique, SEM_oblique,...
                    'LineWidth',3,'Color',bpGlobal.color_list.color_oblique);
        end
        %legend(h,'Oblique task')
end


xlim([0.5,nTimebin + 0.5]);
set(gca, 'fontsize', plotOptions.ftsize)
set(gca, 'TickLabelInterpreter','tex')

ylabel(ylabel_str,'Interpreter','latex')
xlabel(plotOptions.xlabelStr,'Interpreter','latex')
box off
end