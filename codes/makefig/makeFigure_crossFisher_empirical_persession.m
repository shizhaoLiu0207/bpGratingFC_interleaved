clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
%% 1. scatter plot
doThis = 0;
if doThis
    %plotOptions.plotfield = 'fisher';  % 1. I_real v.s. I_cross per session - scatter
    % plotOptions.plotfield = 'deltafisher';  % 2. I redundacy within v.s cross per session - scatter
    plotOptions.plotfield = 'deltafisher_percent'; % 3. I redundacy (percent) within v.s cross per session - scatter
    plotOptions.plottype = 'scatter';
    switch plotOptions.plotfield
        case 'fisher'
            xlabelStr = '$I_\textrm{real}$';
            ylabelStr = '$I_\textrm{cross}$';
        case 'deltafisher'
            xlabelStr = '$I_\textrm{redundancy,real}$';
            ylabelStr = '$I_\textrm{redundancy,cross}$';
        case 'deltafisher_percent'
            
            xlabelStr = '$I_\textrm{redundancy,real}$ (Percent)';
            ylabelStr = '$I_\textrm{redundancy,cross}$ (Percent)';
    end
    plotOptions.ftsize = 14;
    figure;
    set(gcf,'unit','inches','position',[0,0,12,4]);
    
    ax_1 = subplot(1,4,1);
    set(ax_1,'position',get(ax_1,'position')+[-0.06 0.05 0.01 -0.08]);
    session_list_plot = bpGlobal.rolo.session_list.switching;
    plotOptions.task = 'cardinal';
    fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
    xlabel(xlabelStr,'Interpreter','latex');
    ylabel(ylabelStr,'Interpreter','latex');
    title('Cardinal','Color',bpGlobal.color_list.color_cardinal);
    
    
    ax_2 = subplot(1,4,2);
    set(ax_2,'position',get(ax_2,'position')+[-0.05 0.05 0.01 -0.08]);
    session_list_plot = bpGlobal.rolo.session_list.switching;
    plotOptions.task = 'oblique';
    fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
    xlabel(xlabelStr,'Interpreter','latex');
    title('Oblique','Color',bpGlobal.color_list.color_oblique);
    
    
    
    ax_3 = subplot(1,4,3);
    set(ax_3,'position',get(ax_3,'position')+[-0.01 0.05 0.01 -0.08]);
    session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
    plotOptions.task = 'cardinal';
    fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
    xlabel(xlabelStr,'Interpreter','latex');
    ylabel(ylabelStr,'Interpreter','latex');
    title('Cardinal','Color',bpGlobal.color_list.color_cardinal);
    
    ax_4 = subplot(1,4,4);
    set(ax_4,'position',get(ax_4,'position')+[-0.01 0.05 0.01 -0.08]);
    session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
    plotOptions.task = 'oblique';
    fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
    xlabel(xlabelStr,'Interpreter','latex');
    title('Oblique','Color',bpGlobal.color_list.color_oblique);
    
    save_folder = '../../figures/figures_final/fisher_info_session';
    save_name = fullfile(save_folder, sprintf('scatter_session_%s.svg',plotOptions.plotfield));
    print(save_name,'-dsvg','-vector');
end


%%
% %plotOptions.plotfield = 'fisher';  % 1. I_real v.s. I_cross per session - scatter
% % plotOptions.plotfield = 'deltafisher';  % 2. I redundacy within v.s cross per session - scatter
%  plotOptions.plotfield = 'deltafisher_percent'; % 3. I redundacy (percent) within v.s cross per session - scatter
% plotOptions.plottype = 'timecourse';
% % switch plotOptions.plotfield
% %     case 'fisher'
% %         xlabelStr = '$I_\textrm{real}$';
% %         ylabelStr = '$I_\textrm{cross}$';
% %     case 'deltafisher'
% %         xlabelStr = '$I_\textrm{redundancy,real}$';
% %         ylabelStr = '$I_\textrm{redundancy,cross}$';
% %     case 'deltafisher_percent'
% % 
% %         xlabelStr = '$I_\textrm{redundancy,real}$ (Percent)';
% %         ylabelStr = '$I_\textrm{redundancy,cross}$ (Percent)';
% % end
% plotOptions.ftsize = 14;
% figure;
% set(gcf,'unit','inches','position',[0,0,12,4]);
% 
% ax_1 = subplot(1,2,1);
% set(ax_1,'position',get(ax_1,'position')+[-0.06 0.05 0.01 -0.08]);
% session_list_plot = bpGlobal.rolo.session_list.switching;
% %plotOptions.task = 'cardinal';
% fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
% xlabel(xlabelStr,'Interpreter','latex');
% ylabel(ylabelStr,'Interpreter','latex');
% title('Cardinal','Color',bpGlobal.color_list.color_cardinal);
% 
% 
% % ax_2 = subplot(1,4,2);
% % set(ax_2,'position',get(ax_2,'position')+[-0.05 0.05 0.01 -0.08]);
% session_list_plot = bpGlobal.rolo.session_list.switching;
% plotOptions.task = 'oblique';
% fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
% xlabel(xlabelStr,'Interpreter','latex');
% title('Oblique','Color',bpGlobal.color_list.color_oblique);
% 
% 
% 
% ax_2 = subplot(1,2,2);
% set(ax_2,'position',get(ax_2,'position')+[-0.01 0.05 0.01 -0.08]);
% session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
% plotOptions.task = 'cardinal';
% fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
% xlabel(xlabelStr,'Interpreter','latex');
% ylabel(ylabelStr,'Interpreter','latex');
% title('Cardinal','Color',bpGlobal.color_list.color_cardinal);
% 
% % ax_4 = subplot(1,4,4);
% % set(ax_4,'position',get(ax_4,'position')+[-0.01 0.05 0.01 -0.08]);
% session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
% plotOptions.task = 'oblique';
% fig_it.plot_crossInfo_session_scatter(results_all, session_list_plot, plotOptions);
% xlabel(xlabelStr,'Interpreter','latex');
% title('Oblique','Color',bpGlobal.color_list.color_oblique);
% 
% save_folder = '../../figures/figures_final/fisher_info_session';
% save_name = fullfile(save_folder, sprintf('scatter_session_%s.svg',plotOptions.plotfield));
% print(save_name,'-dsvg','-vector');
