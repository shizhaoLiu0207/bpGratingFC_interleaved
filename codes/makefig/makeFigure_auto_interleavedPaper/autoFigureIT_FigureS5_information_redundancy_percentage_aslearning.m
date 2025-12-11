clear all
clc
close all

global   bpGlobal ftsize
ftsize                                     = 14;
bpGratingFCGlobal();


versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';

[session_list, results_neural] = data_prep(versionName);
%%
figure
set(gcf,'unit','inches','position',[0,0,12,8])
save_folder = '../../../figures/figures_auto_interleavedpaper';
save_name = fullfile(save_folder,'Figure_S5_I_redundancy_percentage_interleaved_sizeControl_interleaved.svg');
tex_name = fullfile(save_folder,'Figure_S5_I_redundancy_percentage_interleaved_sizeControl_report.tex');

tex_bar_name = fullfile(save_folder,'Figure_S5_I_redundancy_percentage_interleaved_sizeControl_report_bar.tex');
%% 1. I_redundancy timecourse - Monkey R
ax_1 = subplot(2,4,[1,2]);
set(ax_1,'position',get(ax_1,'position')+[-0.04 0 0.03 0]);
plotOptions.dataName = 'combine_delta_percent';
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.plotInterleaved = 0;
plotOptions.uncertaintyName = 'CI';
plotOptions.plotZeroLine = true;
% plotOptions.plotOneLine = true;

fig_it.plot_timecourse_interleaved_outside(results_neural,'Ro',session_list.rolo_switching,plotOptions)
ylim([-0.3,1] * 100)
set(gca,'TickDir','out','fontsize',ftsize)
ylabel('$I_\textrm{redundancy}$(percentage)','Interpreter','latex','fontsize',ftsize+2)
title('Monkey R', 'fontsize',ftsize+2)
% text(6,0.8*100,'Cardinal Epoch','color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
% text(30,0.8*100,'Oblique Epoch','color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');

%% 2. I_redundancy timecourse - Monkey G
ax_2 = subplot(2,4,[3,4]);
set(ax_2,'position',get(ax_2,'position')+[0.02 0 0.03 0]);
plotOptions.dataName = 'combine_delta_percent';
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.plotInterleaved = 0;
plotOptions.uncertaintyName = 'CI';
fig_it.plot_timecourse_interleaved_outside(results_neural,'Gr',session_list.gremlin_switching,plotOptions)
ylim([-2,1.5] * 100)

% text(60,100,'Cardinal Epoch','color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
% text(15,100,'Oblique Epoch','color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');
set(gca,'TickDir','out','fontsize',ftsize)
title('Monkey G', 'fontsize',ftsize+2)
ylabel('$I_\textrm{redundancy}$(percentage)','Interpreter','latex','fontsize',ftsize+2)

%% 3. Scatter (I_redundancy vs learning index) - Monkey R cardinal
ax_3 = subplot(2,4,5);
set(ax_3,'position',get(ax_3,'position')+[-0.04 0 0.01 -0.04]);
plotOptions.dataStr = 'combine_delta_percent';
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.plotInterleaved = 0;
plotOptions.doPartial = 0;
plotOptions.printStr = 0;

plotOptions.epoch='cardinalInterleaved';%
% 'obliqueOnly'};
plotOptions.subjectCode = 'Ro';
stats_info = fig.plot_behav_fisher_scatterplot(results_neural,plotOptions);
%print(fullfile(saveFolder,'scatter_subset_monkeyR_paper.svg'), '-dsvg');
xlabel('Learning Index')
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
set(gca,'TickDir','out','fontsize',ftsize)
text(0.035,-0.15 * 100, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r, fig.p2star(stats_info.p)),...
    'Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'Interpreter','latex')

ylabel('$I_\textrm{redundancy}$(percentage)','Interpreter','latex','fontsize',ftsize+2)

monkeyR_cardinal_text = sprintf(' Monkey R, cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$', stats_info.df, stats_info.r, stats_info.p);
%% 4. Scatter (I_redundancy vs learning index) - Monkey R oblique
ax_4 = subplot(2,4,6);
set(ax_4,'position',get(ax_4,'position')+[-0.03 0 0.01 -0.04]);
plotOptions.epoch   =   'obliqueInterleaved';
plotOptions.subjectCode = 'Ro';
stats_info = fig.plot_behav_fisher_scatterplot(results_neural,plotOptions);
ylabel(''); xlabel('Learning Index')
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');
set(gca,'TickDir','out','fontsize',ftsize)
text(0.035,50, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r, fig.p2star(stats_info.p)),...
    'Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'Interpreter','latex')

monkeyR_oblique_text = sprintf('Monkey R, oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$', stats_info.df, stats_info.r, stats_info.p);

%% 5. Scatter (I_redundancy vs learning index) - Monkey G cardinal
ax_5 = subplot(2,4,7);
set(ax_5,'position',get(ax_5,'position')+[0.01 0 0.01 -0.04]);
plotOptions.dataStr = 'combine_delta_percent';
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.plotInterleaved = 0;
plotOptions.doPartial = 0;

plotOptions.epoch = 'cardinalInterleaved';
plotOptions.subjectCode = 'Gr';
stats_info = fig.plot_behav_fisher_scatterplot(results_neural,plotOptions);
 xlabel('Learning Index')
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
set(gca,'TickDir','out','fontsize',ftsize)
text(0.06,-1.2*100, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r, fig.p2star(stats_info.p)),...
    'Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'Interpreter','latex')
ylabel('$I_\textrm{redundancy}$(percentage)','Interpreter','latex','fontsize',ftsize+2)

monkeyG_oblique_text = sprintf(' Monkey G, oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$', stats_info.df, stats_info.r, stats_info.p);
%% 6. Scatter (I_redundancy vs learning index) - Monkey G oblique
ax_6 = subplot(2,4,8);
set(ax_6,'position',get(ax_6,'position')+[0.03 0 0.01 -0.04]);
plotOptions.epoch = 'obliqueInterleaved';
plotOptions.subjectCode = 'Gr';
stats_info = fig.plot_behav_fisher_scatterplot(results_neural,plotOptions);
ylabel(' ');
xlabel('Learning Index')
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');
set(gca,'TickDir','out','fontsize',ftsize)
text(0.035,0.80*100, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r, fig.p2star(stats_info.p)),...
    'Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'Interpreter','latex')

monkeyG_cardinal_text = sprintf(' Monkey G, cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$', stats_info.df, stats_info.r, stats_info.p);

%%
annotation('textbox',[0,0.95,0.1,0.04],'string','a','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0,0.45,0.1,0.04],'string','b','fontsize',40,'FontWeight','bold','EdgeColor','none')
%% save figure as svg
print(save_name,'-dsvg','-vector')

%% write statistic report
fid = fopen(tex_name,'wt');

fwrite(fid, ['Spearman rank correlation',...
    monkeyR_cardinal_text, monkeyR_oblique_text, monkeyG_oblique_text, monkeyG_cardinal_text]);

fclose(fid)


%%

function results_neural = data_prep_passive(results_neural)
idx_remove = [];
session_str_list = unique({results_neural(:).sessionStr});
for n = 1:numel(session_str_list)
    idx = strcmp({results_neural(:).sessionStr}, session_str_list{n});
    if ~ismember('mainTask', {results_neural(idx).sessionType})
        idx_remove = [idx_remove, find(idx)];
    end
end
results_neural(idx_remove) = [];


%idx_outlier_session = [337];
idx_outlier_session = find(strcmp({results_neural(:).sessionStr},'Ro20211220'));
for i = 1:numel(idx_outlier_session)
    n = idx_outlier_session(i);
    %X = results_sizeControl_combined(n).combine_delta_sample;
    idx_outlier_delta = detect_outlier(results_neural(n).combine_delta_sample);
    idx_outlier_fisher = detect_outlier(results_neural(n).combine_fisher_sample);
    idx_outlier_shuffle = detect_outlier(results_neural(n).combine_shuffle_sample);

    idx_remove = idx_outlier_delta | idx_outlier_fisher | idx_outlier_shuffle;


    results_neural(n).combine_delta_sample(idx_remove) = [];
    results_neural(n).combine_fisher_sample(idx_remove) = [];
    results_neural(n).combine_shuffle_sample(idx_remove) = [];

    results_neural(n).combine_delta = mean(results_neural(n).combine_delta_sample);
    results_neural(n).combine_delta_std = std(results_neural(n).combine_delta_sample);

    results_neural(n).combine_fisher = mean(results_neural(n).combine_fisher_sample);
    results_neural(n).combine_fisher_std = std(results_neural(n).combine_fisher_sample);

    results_neural(n).combine_shuffle = mean(results_neural(n).combine_shuffle_sample);
    results_neural(n).combine_shuffle_std = std(results_neural(n).combine_shuffle_sample);
end


end

function idx_outlier = detect_outlier(X)

idx_outlier = isoutlier(X,'ThresholdFactor',3) & isoutlier(X,'percentile',[1,99]);

end