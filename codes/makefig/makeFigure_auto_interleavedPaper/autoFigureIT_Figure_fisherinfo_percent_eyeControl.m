clear all
clc
close all

global   bpGlobal ftsize
ftsize = 14;
bpGratingFCGlobal();

session_list_rolo      = bpGlobal.rolo.session_list.switching;
session_list_gremlin   = bpGlobal.gremlin.session_list.interleaved;
%%
eye_var_option_list = {'position';'velocity';'pupil'};
task_list_list = {{'cardinal'};{'oblique'};{'cardinal'};{'oblique'}};
subjectCode_list = {'Ro';'Ro';'Gr';'Gr'};
position_col_shift_list = {[-0.05,0,0.03,0];[-0.04,0,0.03,0];...
    [0.01,0,0.03,0];[0.04,0,0.03,0]};
position_row_shift_list = {[0,-0.02,0,0];[0,-0.03,0,0];[0,-0.04,0,0]};
ylim_list = {[-5,25];[-5,10];[-5,50];[-10,15]};
figure
set(gcf,'unit','inches','position',[0,0,12,10])
for i_row = 1:3
    eye_var_option = eye_var_option_list{i_row};
    for i_col = 1:4
        task_list = task_list_list{i_col};
        subjectCode = subjectCode_list{i_col};

        ax_1 = subplot(3,4,(i_row - 1) * 4 + i_col);
        set(ax_1,'position',get(ax_1,'position')+position_col_shift_list{i_col}+ position_row_shift_list{i_row});
        switch subjectCode
            case 'Ro'
                session_list_plot = session_list_rolo;
            case 'Gr'
                session_list_plot = session_list_gremlin;
        end
        switch eye_var_option
            case 'position'
                filter_name_all  = 'all_trials_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_75_small = 'distanceEyePosition_center_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_75_large = 'distanceEyePosition_surround_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_50_small = 'distanceEyePosition_center_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_50_large = 'distanceEyePosition_surround_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
            case 'velocity'
                filter_name_all  = 'all_trials_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_75_small = 'eyeVelocity_low_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_75_large = 'eyeVelocity_high_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_50_small = 'eyeVelocity_low_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_50_large = 'eyeVelocity_high_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
            case 'pupil'
                filter_name_all  = 'all_trials_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_75_small = 'pupilSize_low_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_75_large = 'pupilSize_high_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_50_small = 'pupilSize_low_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
                filter_name_50_large = 'pupilSize_high_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';
        end


        plotOptions = struct();
        plotOptions.errorbar = 'SEM_session';
        plotOptions.ftsize = 16;
        plotOptions.plotdata = 'delta';
        plotOptions.markersize = 6;
        plotOptions.task_list  = task_list;
        plotOptions.plot_dash_line = false;

        % 1. all trials
        results_all = load_fisher_results(filter_name_all);
        plotOptions.x_cardinal = 1; plotOptions.x_oblique = 1;
        [~,data_all_cardinal, data_all_oblique] = fig_it.plot_diff_errorbar(results_all, session_list_plot, plotOptions);

        % 2. 75%, small
        results_all = load_fisher_results(filter_name_75_small);
        plotOptions.x_cardinal = 1.9;plotOptions.x_oblique = 1.9;
        [~,data_75_small_cardinal, data_75_small_oblique] = fig_it.plot_diff_errorbar(results_all, session_list_plot, plotOptions);

        % 3. 75%, large
        results_all = load_fisher_results(filter_name_75_large);
        plotOptions.x_cardinal = 2.1;plotOptions.x_oblique = 2.1;
        [~,data_75_large_cardinal,data_75_large_oblique,h_c,h_o] = fig_it.plot_diff_errorbar(results_all, session_list_plot, plotOptions);
        if contains(class(h_c),'ErrorBar'), h_c.Bar.LineStyle = 'dotted'; end
        if contains(class(h_o),'ErrorBar'), h_o.Bar.LineStyle = 'dotted'; end

        % 4. 50%, small
        results_all = load_fisher_results(filter_name_50_small);
        plotOptions.x_cardinal = 2.9;plotOptions.x_oblique = 2.9;
        [~,data_50_small_cardinal, data_50_small_oblique] = fig_it.plot_diff_errorbar(results_all, session_list_plot, plotOptions);

        % 5. 50%, large
        results_all = load_fisher_results(filter_name_50_large);
        plotOptions.x_cardinal = 3.1;plotOptions.x_oblique = 3.1;
        [~,data_50_large_cardinal,data_50_large_oblique,h_c,h_o] = fig_it.plot_diff_errorbar(results_all, session_list_plot, plotOptions);
        if contains(class(h_c),'ErrorBar'), h_c.Bar.LineStyle = 'dotted'; end
        if contains(class(h_o),'ErrorBar'), h_o.Bar.LineStyle = 'dotted'; end

        set(gca,'xtick',[1,2,3],'xticklabels',{'All','75%','50%'})
        line([0.8,3.2],[0,0],'linestyle','--','color','black','linewidth',2)

        if ismember('cardinal',task_list)
            offset_coef = 0.2;
            if ttest(data_all_cardinal,data_75_small_cardinal)
                fig.show_ttest(data_all_cardinal, data_75_small_cardinal, [1,1.9], offset_coef);
            end
            if ttest(data_all_cardinal, data_50_small_cardinal)
                fig.show_ttest(data_all_cardinal, data_50_small_cardinal, [1,2.9], offset_coef);
            end
            offset_coef = 0.1;
            if ttest(data_75_small_cardinal, data_75_large_cardinal)
                fig.show_ttest(data_75_small_cardinal, data_75_large_cardinal, [1.9,2.1], offset_coef);
            end
            if ttest(data_50_small_cardinal,data_50_large_cardinal)
                fig.show_ttest(data_50_small_cardinal, data_50_large_cardinal, [2.9,3.1], offset_coef);
            end
            [h_all,p_all,~,stats_all] = ttest(data_all_cardinal);
            [h_75_small,p_75_small,~,stats_75_small] = ttest(data_75_small_cardinal);
            [h_75_large,p_75_large,~,stats_75_large] = ttest(data_75_large_cardinal);
            [h_50_small,p_50_small,~,stats_50_small] = ttest(data_50_small_cardinal);
            [h_50_large,p_50_large,~,stats_50_large] = ttest(data_50_large_cardinal);
        end

        if ismember('oblique',task_list)
            offset_coef = 0.2;
            if ttest(data_all_oblique, data_75_small_oblique)
                fig.show_ttest(data_all_oblique, data_75_small_oblique, [1,1.9], offset_coef);
            end
            if ttest(data_all_oblique, data_50_small_oblique)
                fig.show_ttest(data_all_oblique, data_50_small_oblique, [1,2.9], offset_coef);
            end
            offset_coef = 0.1;
            if ttest(data_75_small_oblique, data_75_large_oblique)
                fig.show_ttest(data_75_small_oblique, data_75_large_oblique, [1.9,2.1], offset_coef);
            end
            if ttest(data_50_small_oblique, data_50_large_oblique)
                fig.show_ttest(data_50_small_oblique, data_50_large_oblique, [2.9,3.1], offset_coef);
            end
            [h_all,p_all,~,stats_all] = ttest(data_all_oblique);
            [h_75_small,p_75_small,~,stats_75_small] = ttest(data_75_small_oblique);
            [h_75_large,p_75_large,~,stats_75_large] = ttest(data_75_large_oblique);
            [h_50_small,p_50_small,~,stats_50_small] = ttest(data_50_small_oblique);
            [h_50_large,p_50_large,~,stats_50_large] = ttest(data_50_large_oblique);
        end

        
        if ismember('oblique',task_list)
            ylabel('')
        end
        ylim(ylim_list{i_col});
        if i_row == 1 & ismember('cardinal',task_list)
            title('Cardinal','color','red')
        end

        if i_row == 1 & ismember('oblique',task_list)
            title('Oblique','color','blue')
        end
    end
end
%%
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0.22,0.95,0.15,0.04],'string','Monkey R','fontsize',20,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.69,0.95,0.15,0.04],'string','Monkey G','fontsize',20,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.01,0.95,0.15,0.04],'string','a','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.01,0.65,0.15,0.04],'string','b','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.01,0.32,0.15,0.04],'string','c','fontsize',40,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.02,0.89,0.15,0.04],'string','Position','fontsize',16,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.02,0.59,0.15,0.04],'string','Velocity','fontsize',16,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.02,0.26,0.15,0.04],'string','Pupil','fontsize',16,'FontWeight','bold','EdgeColor','none')

save_folder = '../../../figures/figures_auto_interleavedpaper';
save_name = fullfile(save_folder,'figure_fisherInfo_cross_percent_eyecontrol.svg');
print(save_name,'-dsvg','-vector')
%% helper functions
function results_all = load_fisher_results(filter_name)


saveFolder = sprintf('../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);
load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);


end