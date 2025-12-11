function  stats_info = plot_info_withintrial(results_all, session_list_plot, plotOptions)

global bpGlobal ftsize
ftsize = plotOptions.ftsize;
nTimebin = 8;
[data_plot_cardinal_within, data_plot_cardinal_cross, data_plot_oblique_within, data_plot_oblique_cross] = deal(cell(nTimebin,1));
for t = 1:nTimebin
    idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == t &...
        ismember({results_all(:).sessionStr}, session_list_plot);

    % 
    % diff_info_percent_cardinal  = cell2mat({results_all(idx).diff_info_percent_cardinal_median});
    % diff_info_percent_oblique   = cell2mat({results_all(idx).diff_info_percent_oblique_median});
    % 
    % diff_delta_percent_cardinal =  cell2mat({results_all(idx).diff_delta_percent_cardinal_median});
    % diff_delta_percent_oblique  =  cell2mat({results_all(idx).diff_delta_percent_oblique_median});

    
    I_delta_cardinal_within         = cell2mat({results_all(idx).delta_cardinal_cardinal_median});
    I_delta_cardinal_cross          = cell2mat({results_all(idx).delta_cardinal_oblique_median});

    I_delta_oblique_within          = cell2mat({results_all(idx).delta_oblique_oblique_median});
    I_delta_oblique_cross           = cell2mat({results_all(idx).delta_oblique_cardinal_median});

    I_delta_percent_cardinal_within =  cell2mat({results_all(idx).delta_percent_cardinal_cardinal_median});
    I_delta_percent_cardinal_cross  =  cell2mat({results_all(idx).delta_percent_cardinal_oblique_median});
    I_delta_percent_oblique_within  =  cell2mat({results_all(idx).delta_percent_oblique_oblique_median});
    I_delta_percent_oblique_cross   =  cell2mat({results_all(idx).delta_percent_oblique_cardinal_median});

    I_cardinal_within               = cell2mat({results_all(idx).fisher_cardinal_cardinal_median});
    I_cardinal_cross                = cell2mat({results_all(idx).fisher_cardinal_oblique_median});
    I_oblique_within                = cell2mat({results_all(idx).fisher_oblique_oblique_median});
    I_oblique_cross                 = cell2mat({results_all(idx).fisher_oblique_cardinal_median});



    switch plotOptions.plotdata
        case 'delta'
            data_plot_cardinal_within{t}    = I_delta_cardinal_within;
            data_plot_cardinal_cross{t}     = I_delta_cardinal_cross;
        
            data_plot_oblique_within{t}     = I_delta_oblique_within;
            data_plot_oblique_cross{t}      = I_delta_oblique_cross;
            ylabel_str = '$I_\textrm{redundancy}$';
        case 'delta_percent'
            data_plot_cardinal_within{t}    = I_delta_percent_cardinal_within;
            data_plot_cardinal_cross{t}     = I_delta_percent_cardinal_cross;
        
            data_plot_oblique_within{t}     = I_delta_percent_oblique_within;
            data_plot_oblique_cross{t}      = I_delta_percent_oblique_cross;
            ylabel_str = '$I_\textrm{redundancy}$(Percent)';
        case 'info'
            data_plot_cardinal_within{t}    = I_cardinal_within;
            data_plot_cardinal_cross{t}     = I_cardinal_cross;
        
            data_plot_oblique_within{t}     = I_oblique_within;
            data_plot_oblique_cross{t}      = I_oblique_cross;
            ylabel_str = 'Linear Fisher information';
    end

    % if strcmp(plotOptions.errorbar,'CI_sample')
    %     I_delta_within_cardinal_CI  = results_all(idx).diff_info_percent_cardinal_CI;
    %     diff_info_percent_oblique_CI  = results_all(idx).diff_info_percent_oblique_CI;
    % 
    %     diff_delta_percent_cardinal_CI =  results_all(idx).diff_delta_percent_cardinal_CI;
    %     diff_delta_percent_oblique_CI  =  results_all(idx).diff_delta_percent_oblique_CI;
    % 
    % 
    % end
  
end

 xlabelStrs = cell(nTimebin,2);
for i = 1:nTimebin
    idx = find(cell2mat({results_all(:).timeWinIndex}) == i, 1);

        %xlabelStrs{i} = num2str(results_Struct(idx).timeWin(2));
        xlabelStrs{i,1} = num2str(results_all(idx).timeWin(1) - 1);
        xlabelStrs{i,2} = num2str(results_all(idx).timeWin(2));
   
end

nSession_cardinal       = cellfun(@numel, data_plot_cardinal_within);
nSession_oblique        = cellfun(@numel, data_plot_oblique_within);
mean_cardinal_within           = cellfun(@mean, data_plot_cardinal_within);
mean_cardinal_cross           = cellfun(@mean, data_plot_cardinal_cross);
mean_oblique_within            = cellfun(@mean, data_plot_oblique_within);
mean_oblique_cross            = cellfun(@mean, data_plot_oblique_cross);
SEM_cardinal_within            = cellfun(@std, data_plot_cardinal_within) ./ sqrt(nSession_cardinal);
SEM_cardinal_cross            = cellfun(@std, data_plot_cardinal_cross) ./ sqrt(nSession_cardinal);
SEM_oblique_within             = cellfun(@std, data_plot_oblique_within) ./ sqrt(nSession_oblique);
SEM_oblique_cross             = cellfun(@std, data_plot_oblique_cross) ./ sqrt(nSession_oblique);

if strcmp(plotOptions.errorbar,'CI_sample')
    CI_cardinal             = cat(1,CI_plot_cardinal{:});
    CI_oblique              = cat(1,CI_plot_oblique{:});
end


%line([0.5,nTimebin + 0.5], [0,0],'color','black','linestyle','--'); hold on
switch plotOptions.task
    case 'cardinal'
        switch plotOptions.errorbar
            case 'CI_sample'
               h = errorbar([1:nTimebin], mean_cardinal, mean_cardinal - CI_cardinal(:,1),  CI_cardinal(:,2) - mean_cardinal,...
                    'LineWidth',3,'Color',bpGlobal.color_list.color_cardinal);
            case 'SEM_session'
                % h(1) = errorbar([1:nTimebin], mean_cardinal_within, SEM_cardinal_within,...
                %     'LineWidth',3,'Color',bpGlobal.color_list.color_cardinal,'LineStyle','-');
                % 
                % h(2) = errorbar([1:nTimebin], mean_cardinal_cross, SEM_cardinal_cross,...
                %     'LineWidth',3,'Color',bpGlobal.color_list.color_cardinal,'LineStyle','--');
                plotMarker = '^';
                %xticklabelStr = {'50','250','450','650','850','1250','1450','1650'};
                plotStyle = '-';
                includeOri = false;
                [h, beta_within, beta_p_within] = fig.plot_errobar_ttest(data_plot_cardinal_within, ...
                    bpGlobal.color_list.color_cardinal,plotMarker, xlabelStrs,ylabel_str, plotStyle, includeOri);

                plotStyle = '--';
                [h, beta_cross, beta_p_cross] = fig.plot_errobar_ttest(data_plot_cardinal_cross, ...
                    bpGlobal.color_list.color_cardinal,plotMarker, xlabelStrs,ylabel_str, plotStyle, includeOri);
                lgd = findobj('type', 'legend');
                set(lgd, 'visible', 'off')
        end
        %egend(h,'Cardinal task')
    case 'oblique'
        % switch plotOptions.errorbar
        %     case 'CI_sample'
        %         h = errorbar([1:nTimebin], mean_oblique, mean_oblique - CI_oblique(:,1),  CI_oblique(:,2) - mean_oblique,...
        %             'LineWidth',3,'Color',bpGlobal.color_list.color_oblique);
        %     case 'SEM_session'
        %         h(1) = errorbar([1:nTimebin], mean_oblique_within, SEM_oblique_within,...
        %             'LineWidth',3,'Color',bpGlobal.color_list.color_oblique,'LineStyle','-');
        %         h(2) = errorbar([1:nTimebin], mean_oblique_cross, SEM_oblique_cross,...
        %             'LineWidth',3,'Color',bpGlobal.color_list.color_oblique,'LineStyle','--');
        % 
        % end
        %legend(h,'Oblique task')
end

stats_info.beta_within = beta_within;
stats_info.beta_p_within = beta_p_within;
stats_info.beta_cross = beta_cross;
stats_info.beta_p_cross = beta_p_cross;
% xlim([0.5,nTimebin + 0.5]);
set(gca, 'fontsize', plotOptions.ftsize)
set(gca, 'TickLabelInterpreter','tex')
% legend(h,'Within','Cross','Location','northwest')
ylabel(ylabel_str,'Interpreter','latex')
xlabel(plotOptions.xlabelStr,'Interpreter','latex')
% xlabel(plotOptions.xlabelStr,'Interpreter','latex')
% box off
% set(gca,'xtick',[0.5:1:8.5],'XTicklabels',{'50','250','450','650','850','1250','1450','1650'})
end