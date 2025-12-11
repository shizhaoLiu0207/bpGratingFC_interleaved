function plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions)
global bpGlobal
 idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
        ismember({results_all(:).sessionStr}, session_list_plot);



nSession = sum(idx);

switch plotOptions.plotfield
    case 'fisher'
        I_cardinal          = cell2mat({results_all(idx).fisher_cardinal_cardinal_median})';
        I_cardinal_cross    = cell2mat({results_all(idx).fisher_cardinal_oblique_median})';
        I_oblique           = cell2mat({results_all(idx).fisher_oblique_oblique_median})';
        I_oblique_cross     = cell2mat({results_all(idx).fisher_oblique_cardinal_median})';

        tmp = {results_all(idx).fisher_cardinal_cardinal_CI};
        I_cardinal_CI          = cat(1,tmp{:});
        tmp = {results_all(idx).fisher_cardinal_oblique_CI};
        I_cardinal_cross_CI    = cat(1,tmp{:});
        tmp = {results_all(idx).fisher_oblique_oblique_CI};
        I_oblique_CI           = cat(1,tmp{:});
        tmp = {results_all(idx).fisher_oblique_cardinal_CI};
        I_oblique_cross_CI     = cat(1,tmp{:});

    case 'deltafisher'
        I_cardinal          = cell2mat({results_all(idx).delta_cardinal_cardinal_median})';
        I_cardinal_cross    = cell2mat({results_all(idx).delta_cardinal_oblique_median})';
        I_oblique           = cell2mat({results_all(idx).delta_oblique_oblique_median})';
        I_oblique_cross     = cell2mat({results_all(idx).delta_oblique_cardinal_median})';
        % 
        tmp = {results_all(idx).delta_cardinal_cardinal_CI};
        I_cardinal_CI          = cat(1,tmp{:});
        tmp = {results_all(idx).delta_cardinal_oblique_CI};
        I_cardinal_cross_CI    = cat(1,tmp{:});
        tmp = {results_all(idx).delta_oblique_oblique_CI};
        I_oblique_CI           = cat(1,tmp{:});
        tmp = {results_all(idx).delta_oblique_cardinal_CI};
        I_oblique_cross_CI     = cat(1,tmp{:});
    case 'deltafisher_percent'
        I_cardinal          = cell2mat({results_all(idx).delta_percent_cardinal_cardinal_median})';
        I_cardinal_cross    = cell2mat({results_all(idx).delta_percent_cardinal_oblique_median})';
        I_oblique           = cell2mat({results_all(idx).delta_percent_oblique_oblique_median})';
        I_oblique_cross     = cell2mat({results_all(idx).delta_percent_oblique_cardinal_median})';

        tmp = {results_all(idx).delta_percent_cardinal_cardinal_CI};
        I_cardinal_CI          = cat(1,tmp{:});
        tmp = {results_all(idx).delta_percent_cardinal_oblique_CI};
        I_cardinal_cross_CI    = cat(1,tmp{:});
        tmp = {results_all(idx).delta_percent_oblique_oblique_CI};
        I_oblique_CI           = cat(1,tmp{:});
        tmp = {results_all(idx).delta_percent_oblique_cardinal_CI};
        I_oblique_cross_CI     = cat(1,tmp{:});
end

%figure
switch plotOptions.plottype 
    case  'scatter'
        switch plotOptions.task
            case 'cardinal'
                   % scatter(I_cardinal, I_cardinal_cross)
                    plot_color          = bpGlobal.color_list.color_cardinal;
                    plot_color_light    = bpGlobal.color_list.color_cardinal_light;
                  
                    plot(I_cardinal, I_cardinal_cross,'.','color',plot_color,'markersize',20); 
                    %plot_scatter_errorbar_cross(I_cardinal,I_cardinal_cross, plot_color, plot_color_light);
                    [~,p,~,stats] = ttest(I_cardinal, I_cardinal_cross);
                    text(max(I_cardinal) * 0.5, min(I_cardinal_cross), sprintf('$t(%d) = %.2f^{%s}$',stats.df, stats.tstat, fig.p2star(p)),...
                        'Interpreter','latex','FontSize',plotOptions.ftsize ,'color',plot_color)
                    lower_lim = min([I_cardinal;I_cardinal_cross],[],'all');
                    upper_lim = max([I_cardinal;I_cardinal_cross],[],'all');
                    line([lower_lim,upper_lim],[lower_lim,upper_lim],'linewidth',3,'color','black','linestyle','--')
            case 'oblique'
                    plot_color          = bpGlobal.color_list.color_oblique;
                    plot_color_light    = bpGlobal.color_list.color_oblique_light;
                 
                    [~,p,~,stats] = ttest(I_oblique, I_oblique_cross);
                    plot(I_oblique, I_oblique_cross,'.','color',plot_color,'markersize',20); 
                    text(max(I_cardinal) * 0.5, min(I_cardinal_cross), sprintf('$t(%d) = %.2f^{%s}$',stats.df, stats.tstat, fig.p2star(p)),...
                        'Interpreter','latex','FontSize',plotOptions.ftsize ,'color',plot_color)
                    %plot_scatter_errorbar_cross(I_oblique,I_oblique_cross,plot_color, plot_color_light);
        
        
                    lower_lim = min([I_oblique;I_oblique_cross],[],'all');
                    upper_lim = max([I_oblique;I_oblique_cross],[],'all');
                    line([lower_lim,upper_lim],[lower_lim,upper_lim],'linewidth',3,'color','black','linestyle','--')
        end
    case 'timecourse'

        h(1) = errorbar([1:nSession],I_cardinal, I_cardinal - I_cardinal_CI(:,1), I_cardinal_CI(:,2) - I_cardinal,...
            'color',bpGlobal.color_list.color_cardinal); hold on
        h(2) = errorbar([1:nSession],I_cardinal_cross, I_cardinal_cross - I_cardinal_cross_CI(:,1), I_cardinal_cross_CI(:,2) - I_cardinal_cross,...
            'color',bpGlobal.color_list.color_cardinal_light); hold on

        h(3) = errorbar([1:nSession],I_oblique, I_oblique - I_oblique_CI(:,1), I_oblique_CI(:,2) - I_oblique,...
            'color',bpGlobal.color_list.color_oblique); hold on
        h(4) = errorbar([1:nSession],I_oblique_cross, I_oblique_cross - I_oblique_cross_CI(:,1), I_oblique_cross_CI(:,2) - I_oblique_cross,...
            'color',bpGlobal.color_list.color_oblique_light); 

        legend(h,'Cardinal','Cardinal cross','Oblique','Oblique cross')


end

box off
set(gca,'fontsize',plotOptions.ftsize);

    function plot_scatter_errorbar_cross(x_median,y_median, plot_color,plot_color_light)
    for n = 1:size(x_median,1)

        if x_CI(n,2) < y_CI(n,1) | y_CI(n,2) < x_CI(n,1)
            color_use = plot_color;
        else
            color_use = plot_color_light;
        end
       
        plot(x_median(n),y_median(n),'.','color',plot_color,'markersize',20); hold on
        % eb(1) = errorbar(x_median(n),y_median(n),x_median(n) - x_CI(n,1),x_CI(n,2) - x_median(n), ...
        %     'horizontal','LineStyle','none');
        % eb(2) = errorbar(x_median(n),y_median(n),y_median(n) - y_CI(n,1),y_CI(n,2) - y_median(n), ...
        %     'vertical','LineStyle','none');
        % set(eb, 'color',color_use,'linewidth',2)
        % 


    end
end

% subplot(2,1,1)
% y_mean = deltaI_cardinal;
% y_errorbar_low = y_mean - deltaI_cardinal_CI(:,1);
% y_errorbar_high = deltaI_cardinal_CI(:,2) - y_mean;
% h(1) = errorbar([1:nTimebin],y_mean, y_errorbar_low,y_errorbar_high, ...
%     'LineWidth',2,'color',bpGlobal.color_list.color_cardinal, 'linestyle','-');
% 
% y_mean = deltaI_cardinal_cross;
% y_errorbar_low = y_mean - deltaI_cardinal_cross_CI(:,1);
% y_errorbar_high = deltaI_cardinal_cross_CI(:,2) - y_mean;
% h(2) = errorbar([1:nTimebin],y_mean, y_errorbar_low,y_errorbar_high, ...
%     'LineWidth',2,'color',bpGlobal.color_list.color_cardinal, 'linestyle','--');
% 
% ylabel('I_{redundancy}')
% xlabel('Session number')
% set(gca,'Fontsize',16)
% legend(h,'Within','Cross')
% 
% subplot(2,1,2)
% y_mean = deltaI_oblique;
% y_errorbar_low = y_mean - deltaI_oblique_CI(:,1);
% y_errorbar_high = deltaI_oblique_CI(:,2) - y_mean;
% h(1) = errorbar([1:nTimebin],y_mean, y_errorbar_low,y_errorbar_high, ...
%     'LineWidth',2,'color',bpGlobal.color_list.color_cardinal, 'linestyle','-');
% 
% y_mean = deltaI_oblique_cross;
% y_errorbar_low = y_mean - deltaI_oblique_cross_CI(:,1);
% y_errorbar_high = deltaI_oblique_cross_CI(:,2) - y_mean;
% h(2) = errorbar([1:nTimebin],y_mean, y_errorbar_low,y_errorbar_high, ...
%     'LineWidth',2,'color',bpGlobal.color_list.color_cardinal, 'linestyle','--');
% ylabel('I_{redundancy}')
% xlabel('Session number')
% set(gca,'Fontsize',16)
% legend(h,'Within','Cross')
end