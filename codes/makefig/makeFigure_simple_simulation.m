clear all
clc
close all
%%%%% Make figure of results of a simplied model simulation
global  bpGlobal
bpGratingFCGlobal();
color_C = bpGlobal.color_list.color_cardinal; color_O = bpGlobal.color_list.color_oblique;

%%
b_PF            = 0;
cardinal_delta  = 0.08;
oblique_delta   = 0.08;
cardinal_prior  = 1;
oblique_prior   = 1;
session_name_str_afterlearning = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);  

cardinal_delta  = 0;
oblique_delta   = 0;

session_name_str_beforelearning = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);  

data_name_afterlearning =['synthData_use_interleaved_',session_name_str_afterlearning];
session_name_afterlearning = ['Model_',session_name_str_afterlearning];

data_name_beforelearning =['synthData_use_interleaved_',session_name_str_beforelearning];
session_name_beforelearning = ['Model_',session_name_str_beforelearning];

data_folder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved/batch_1';
fisher_folder = '../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/subset_32_256_random_1000/individual_sessions_cross';

data_afterlearning = load(fullfile(data_folder, data_name_afterlearning));
data_beforelearning = load(fullfile(data_folder, data_name_beforelearning));

fisher_afterlearning = load(fullfile(fisher_folder, session_name_afterlearning));
fisher_beforelearning = load(fullfile(fisher_folder, session_name_beforelearning));

% results_fisher_afterlearning = organize_sample_fisher(fisher_afterlearning.dat_fisher_cross);
% results_fisher_beforelearning = organize_sample_fisher(fisher_beforelearning.dat_fisher_cross);

results_fisher_afterlearning = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_afterlearning.dat_fisher_cross);
results_fisher_beforelearning = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_beforelearning.dat_fisher_cross);


results_fisher_afterlearning = get_sample_CI_cross(results_fisher_afterlearning);
results_fisher_beforelearning = get_sample_CI_cross(results_fisher_beforelearning);
%% 1. psychometric curve
figure
set(gcf,'Units','inches','Position',[0,0,6,4])
plotOptions.style_cardinal = '-';
plotOptions.style_oblique = '-';
plotOptions.ftsize = 14;
fig_it.plot_synth_interleaved_psycurve(data_afterlearning.synthData_interleaved, plotOptions); 

plotOptions.style_cardinal = '--';
plotOptions.style_oblique = '--';
plotOptions.ftsize = 14;
fig_it.plot_synth_interleaved_psycurve(data_beforelearning.synthData_interleaved, plotOptions); 

text(5, 0.60,'Before learning','fontsize',14);
text(5, 0.85,'After learning','fontsize',14);
save_folder = '../../figures/figures_final/model_behav';
save_name = fullfile(save_folder,'model_psycurves_interleaved.svg');
print(save_name,'-dsvg','-vector')

%% 2. compare I_real and I_cross
figure
set(gcf,'unit','inches','position',[0,0,6,4])
ax_1 = subplot(3,2,[1,3]); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 14;
fig_it.plot_bar_cross_Info(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions); 
ylim([0, 0.2])
text(0.5,0.18,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,0.18,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Before learning')

ax_2 = subplot(3,2,5 ); hold on
set(ax_2,'position',get(ax_2,'position')+[-0.02 0 0.03 -0.02]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'info';
plotOptions.markersize = 6;
fig_it.plot_diff_errorbar(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
ylabel('Diff.(\%)','Interpreter','latex')
%title('Before learning')
ylim([-10,55])


ax_3 = subplot(3,2,[2,4]); hold on
set(ax_3,'position',get(ax_3,'position')+[0 0.04 0.02 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 14;
fig_it.plot_bar_cross_Info(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions); 
ylabel('')
title('After learning')
ylim([0, 0.2])
text(0.5,0.18,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,0.18,'Oblique','color','blue','FontSize',14,'FontWeight','bold')

ax_4 = subplot(3,2,6); hold on
set(ax_4,'position',get(ax_4,'position')+[0 0 0.02 -0.02]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'info';
plotOptions.markersize = 6;
fig_it.plot_diff_errorbar(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
ylabel('')
%title('Before learning')
ylim([-10,55])


save_folder = '../../figures/figures_final/model_fisher_ideal';
save_name = fullfile(save_folder,'model_fisher_real_cross.svg');
print(save_name,'-dsvg','-vector')
%% 3. compare I_redundancy_percent and I_redundancy_cross_percent
figure
set(gcf,'unit','inches','position',[0,0,6,4])
ax_1 = subplot(3,2,[1,3]); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 14;
plotOptions.plotPercent = true;
fig_it.plot_bar_cross_deltaInfo(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions); 
ylim([40, 100])
text(0.5,90,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,90,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Before learning')



ax_2 = subplot(3,2,5 ); hold on
set(ax_2,'position',get(ax_2,'position')+[-0.02 0 0.03 -0.02]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
fig_it.plot_diff_errorbar(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
ylabel('Diff.(\%)','Interpreter','latex')
%title('Before learning')
ylim([-5,20])


ax_3 = subplot(3,2,[2,4]); hold on
set(ax_3,'position',get(ax_3,'position')+[0 0.04 0.02 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 14;
plotOptions.plotPercent = true;
fig_it.plot_bar_cross_deltaInfo(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions); 
ylabel('')
ylim([40, 100])
text(0.5,90,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,90,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Before learning')

ax_4 = subplot(3,2,6); hold on
set(ax_4,'position',get(ax_4,'position')+[0 0 0.02 -0.02]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
fig_it.plot_diff_errorbar(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
ylabel('')
%title('Before learning')
ylim([-5,20])

save_folder = '../../figures/figures_final/model_fisher_ideal';
save_name = fullfile(save_folder,'model_redundancy_percent_real_cross.svg');
print(save_name,'-dsvg','-vector')
%% 4. compare I_redundancy and I_redundancy_cross
figure
set(gcf,'unit','inches','position',[0,0,10,4])
ax_1 = subplot(1,2,1); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.01 0.03 0.04 -0.03]);

plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 14;
plotOptions.plotPercent = false;
fig_it.plot_bar_cross_deltaInfo(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions); 

ylim([0, 0.42])
text(0.5,0.38,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,0.38,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Before learning')

ax_2 = subplot(1,2,2); hold on
set(ax_2,'position',get(ax_2,'position')+[0.01 0.03 0.04 -0.03]);
fig_it.plot_bar_cross_deltaInfo(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions); 
ylabel('')
ylim([0, 0.42])
text(0.5,0.38,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,0.38,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('After learning')


save_folder = '../../figures/figures_final/model_fisher_ideal';
save_name = fullfile(save_folder,'model_redundancy_real_cross.svg');
print(save_name,'-dsvg','-vector')
% %% 5. (I_cross - I_real / I_cross)
% figure
% set(gcf,'unit','inches','position',[0,0,4,3])
% ax_1 = subplot(1,2,1); hold on
% set(ax_1,'position',get(ax_1,'position')+[0.01 0.03 0.04 -0.03]);
% plotOptions = struct();
% plotOptions.errorbar = 'CI_sample';
% plotOptions.ftsize = 14;
% plotOptions.plotdata = 'info';
% fig_it.plot_diff_errorbar(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
% title('Before learning')
% ylim([-10,55])
% 
% ax_2 = subplot(1,2,2); hold on
% set(ax_2,'position',get(ax_2,'position')+[-0.01 0.03 0.04 -0.03]);
% fig_it.plot_diff_errorbar(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
% ylabel('')
% set(gca,'YColor', 'none')
% title('After learning')
% ylim([-10,55])
% 
% save_folder = '../../figures/figures_final/model_fisher_ideal';
% save_name = fullfile(save_folder,'model_cross_minus_real_percent.svg');
% print(save_name,'-dsvg','-vector')
% %% 6. I_redundancy_within - I_redundancy_cross
% figure
% set(gcf,'unit','inches','position',[0,0,4,3])
% ax_1 = subplot(1,2,1); hold on
% set(ax_1,'position',get(ax_1,'position')+[0.01 0.03 0.04 -0.03]);
% plotOptions = struct();
% plotOptions.errorbar = 'CI_sample';
% plotOptions.ftsize = 14;
% plotOptions.plotdata = 'delta';
% fig_it.plot_diff_errorbar(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
% title('Before learning')
% ylim([-5,20])
% 
% ax_2 = subplot(1,2,2); hold on
% set(ax_2,'position',get(ax_2,'position')+[-0.01 0.03 0.04 -0.03]);
% fig_it.plot_diff_errorbar(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
% ylabel('')
% set(gca,'YColor', 'none')
% ylim([-5,20])
% title('After learning')
% 
% save_folder = '../../figures/figures_final/model_fisher_ideal';
% save_name = fullfile(save_folder,'diff_redundancy_percent.svg');
% print(save_name,'-dsvg','-vector')
%% helper function

function results_cross_sizeControl_perCohr = organize_sample_fisher(dat_fisher_cross)
cohr_list = unique(cell2mat({dat_fisher_cross(:).coherence_level}));
t = 1;
for i = 1:numel(cohr_list)
    idx =  cell2mat({dat_fisher_cross(:).coherence_level}) == cohr_list(i);
    dat_fisher_cross_sub = dat_fisher_cross(idx);
    
    
    results_cross_sizeControl_perCohr(t).sessionStr     = dat_fisher_cross_sub(1).sessionStr;
    results_cross_sizeControl_perCohr(t).sessionType    = dat_fisher_cross_sub(1).sessionType;
    results_cross_sizeControl_perCohr(t).timeWin        = dat_fisher_cross_sub(1).timeWin;
    results_cross_sizeControl_perCohr(t).timeWinIndex   = dat_fisher_cross_sub(1).timeWinIndex;
    results_cross_sizeControl_perCohr(t).nUnit          = dat_fisher_cross_sub(1).N;
    results_cross_sizeControl_perCohr(t).coherence_level = dat_fisher_cross_sub(1).coherence_level;
    
    
    
    results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_cardinal_bc});
    results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_oblique_bc});
    
    results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_cardinal_bc});
    results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_oblique_bc});
    
    results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_cardinal_shuffle_bc});
    results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_oblique_shuffle_bc});
    
    results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_cardinal_shuffle_bc});
    results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_oblique_shuffle_bc});
    t = t+1;
end
    

end