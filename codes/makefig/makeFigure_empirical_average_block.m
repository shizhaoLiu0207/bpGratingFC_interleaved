clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
%filter_name = 'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl';
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
load(fullfile(saveFolder, 'results_SubsampleOrganized_perCohr_fisherInfo_all_sessions'));
results_all_perCohr = results_cross_sizeControl_perCohr;
results_all_perCohr = get_sample_CI_cross(results_all_perCohr);

%% plot average neural signature of switching across two tasks as a function of block size
rolo_block_size = bpGlobal.rolo.session_list.mini_block;
gremlin_block_size = bpGlobal.gremlin.session_list.interleaved_blockSize;

rolo_mini_block_list = rolo_block_size(cell2mat(rolo_block_size(:,2)) < 200, 1);

rolo_random_list = setdiff(bpGlobal.rolo.session_list.switching, rolo_block_size(:,1));

gremlin_block_list = gremlin_block_size(cell2mat(gremlin_block_size(:,2)) > 0,1);
gremlin_random_list = setdiff(bpGlobal.gremlin.session_list.interleaved_real, gremlin_block_list);


%%% block (G), block (R), mini-block(R), trial-by-trial (G,R)
figure
set(gcf,'unit','inches','position',[0,0,10,8])
for k = 1:2
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 16;
if k == 1
    plotOptions.plotdata = 'info';
else
    plotOptions.plotdata = 'delta';
end
plotOptions.markersize = 10;

subplot(2,1,k)

plotOptions.x_val = -0.2;
plotOptions.subjectCode = 'Ro';
session_list_plot = bpGlobal.rolo.session_list.switching;
[stats_info,h(1)] = fig_it.plot_diff_errorbar_avg(results_all, session_list_plot, plotOptions);

plotOptions.x_val = 0.2;
plotOptions.subjectCode = 'Gr';
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info,h(2)] = fig_it.plot_diff_errorbar_avg(results_all, session_list_plot, plotOptions);

plotOptions.x_val = 1;
plotOptions.subjectCode = 'Gr';
session_list_plot = gremlin_block_list;
stats_info = fig_it.plot_diff_errorbar_avg(results_all, session_list_plot, plotOptions);


plotOptions.x_val = 2;
plotOptions.subjectCode = 'Ro';
session_list_plot = rolo_mini_block_list;
stats_info = fig_it.plot_diff_errorbar_avg(results_all, session_list_plot, plotOptions);

plotOptions.x_val = 2.8;
plotOptions.subjectCode = 'Ro';
session_list_plot = rolo_random_list;
stats_info = fig_it.plot_diff_errorbar_avg(results_all, session_list_plot, plotOptions);

plotOptions.x_val = 3.2;
plotOptions.subjectCode = 'Gr';
session_list_plot = gremlin_random_list;
stats_info = fig_it.plot_diff_errorbar_avg(results_all, session_list_plot, plotOptions);

set(gca, 'xtick',[0,1,2,3],'xticklabels',{'All sessions';'Blockwise';'Mini-block';'Trial-by-trial'})
line([-0.7,3.7],[0,0],'linestyle','--','color','black','linewidth',2)

%legend(h,'Monkey R','Monkey G','Location','northoutside','Orientation','horizontal')

end

save_folder = '../../figures/figures_final/fisher_info_bar';
save_name = fullfile(save_folder,'average_neural_signature_block_trial.svg')
print(save_name,'-dsvg','-vector');