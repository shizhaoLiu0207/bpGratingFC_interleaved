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

data_folder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved';
fisher_folder = '../../results/neural/fisherInfo_direct/fisherInfo_direct_modelInterleaved_versionControl/subset_32_256_random_1000/individual_sessions_cross';

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
%% psychometric curve
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
% %% I_redundancy
% figure; 
% set(gcf,'Units','inches','Position',[0,0,12,4])
% 
% ax_1 = subplot(1,3,1); hold on
% set(ax_1,'position',get(ax_1,'position')+[-0.06 0.03 0.01 -0.03]);
% deltaI_cardinal  = results_fisher_beforelearning.delta_cardinal_cardinal_median;
% deltaI_oblique   = results_fisher_beforelearning.delta_oblique_oblique_median;
% deltaI_cardinal_CI  = results_fisher_beforelearning.delta_cardinal_cardinal_CI;
% deltaI_oblique_CI   = results_fisher_beforelearning.delta_oblique_oblique_CI;
% 
% 
% fig_it.plot_bar_errorbar(1, deltaI_cardinal, deltaI_cardinal_CI, color_C, color_C)
% fig_it.plot_bar_errorbar(2, deltaI_oblique, deltaI_oblique_CI, color_O, color_O)
% 
% deltaI_cardinal  = results_fisher_afterlearning.delta_cardinal_cardinal_median;
% deltaI_cardinal_CI  = results_fisher_afterlearning.delta_cardinal_cardinal_CI;
% deltaI_oblique   = results_fisher_afterlearning.delta_oblique_oblique_median;
% deltaI_oblique_CI   = results_fisher_afterlearning.delta_oblique_oblique_CI;
% 
% fig_it.plot_bar_errorbar(4.5, deltaI_cardinal, deltaI_cardinal_CI, color_C, color_C)
% fig_it.plot_bar_errorbar(5.5, deltaI_oblique, deltaI_oblique_CI, color_O, color_O)
% 
% text(0, 0.4,'Before learning','fontsize',14)
% text(3.5, 0.4,'After learning','fontsize',14)
% set(gca,'xtick',[]); set(gca,'fontsize',14)
% set(gca,'ylim',[0,0.45])
% set(gca,'xtick',[1,2,4,5],'xticklabels',{'\color{red}{Cardinal}','\color{blue}{Oblique}','\color{red}{Cardinal}','\color{blue}{Oblique}'})
% set(gca, 'TickLabelInterpreter','tex')
% 
% ylabel('$I_\textrm{redundancy}$','Interpreter','latex')

%% compare I_real and I_cross
figure
set(gcf,'unit','inches','position',[0,0,6,3])
ax_1 = subplot(1,2,1); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.01 0.03 0.04 -0.03]);

I_real_cardinal = results_fisher_beforelearning.fisher_cardinal_cardinal_median;
I_real_cardinal_CI = results_fisher_beforelearning.fisher_cardinal_cardinal_CI;

I_cross_cardinal = results_fisher_beforelearning.fisher_cardinal_oblique_median;
I_cross_cardinal_CI = results_fisher_beforelearning.fisher_cardinal_oblique_CI;

I_real_oblique = results_fisher_beforelearning.fisher_oblique_oblique_median;
I_real_oblique_CI = results_fisher_beforelearning.fisher_oblique_oblique_CI;

I_cross_oblique = results_fisher_beforelearning.fisher_oblique_cardinal_median;
I_cross_oblique_CI = results_fisher_beforelearning.fisher_oblique_cardinal_CI;

fig_it.plot_bar_errorbar(1, I_real_cardinal, I_real_cardinal_CI, color_C, color_C)
fig_it.plot_bar_errorbar(3, I_cross_cardinal, I_cross_cardinal_CI, color_C, color_O)
fig_it.plot_bar_errorbar(6, I_real_oblique, I_real_oblique_CI, color_O, color_O)
fig_it.plot_bar_errorbar(8, I_cross_oblique, I_cross_oblique_CI, color_O, color_C)

ylim([0, 0.2])
text(0.5,0.18,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,0.18,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
set(gca,'xtick',[]); set(gca,'fontsize',14)
set(gca,'xtick',[1,3,6,8],'xticklabels',{'$I_\textrm{real}$';'$I_\textrm{cross}$';'$I_\textrm{real}$';'$I_\textrm{cross}$'})
ylabel('Linear Fisher information','Interpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
title('Before learning')

ax_2 = subplot(1,2,2); hold on
set(ax_2,'position',get(ax_2,'position')+[0.01 0.03 0.04 -0.03]);
I_real_cardinal = results_fisher_afterlearning.fisher_cardinal_cardinal_median;
I_real_cardinal_CI = results_fisher_afterlearning.fisher_cardinal_cardinal_CI;

I_cross_cardinal = results_fisher_afterlearning.fisher_cardinal_oblique_median;
I_cross_cardinal_CI = results_fisher_afterlearning.fisher_cardinal_oblique_CI;

I_real_oblique = results_fisher_afterlearning.fisher_oblique_oblique_median;
I_real_oblique_CI = results_fisher_afterlearning.fisher_oblique_oblique_CI;

I_cross_oblique = results_fisher_afterlearning.fisher_oblique_cardinal_median;
I_cross_oblique_CI = results_fisher_afterlearning.fisher_oblique_cardinal_CI;


fig_it.plot_bar_errorbar(1, I_real_cardinal, I_real_cardinal_CI, color_C, color_C)
fig_it.plot_bar_errorbar(3, I_cross_cardinal, I_cross_cardinal_CI, color_C, color_O)
fig_it.plot_bar_errorbar(6, I_real_oblique, I_real_oblique_CI, color_O, color_O)
fig_it.plot_bar_errorbar(8, I_cross_oblique, I_cross_oblique_CI, color_O, color_C)

title('After learning')

ylim([0, 0.2])
text(0.5,0.18,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,0.18,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
set(gca,'xtick',[]); set(gca,'fontsize',14)
set(gca,'xtick',[1,3,6,8],'xticklabels',{'$I_\textrm{real}$';'$I_\textrm{cross}$';'$I_\textrm{real}$';'$I_\textrm{cross}$'})
%ylabel('Linear Fisher information','Interpreter','latex')
set(gca, 'TickLabelInterpreter','latex')

save_folder = '../../figures/figures_final/model_fisher';
save_name = fullfile(save_folder,'model_fisher_real_cross.svg');
print(save_name,'-dsvg','-vector')
%%
%%% compare I_redundancy and I_redundancy_cross
figure
set(gcf,'unit','inches','position',[0,0,6,3])
ax_1 = subplot(1,2,1); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.01 0.03 0.04 -0.03]);

deltaI_cardinal             = results_fisher_beforelearning.delta_cardinal_cardinal_median;
deltaI_cardinal_CI          = results_fisher_beforelearning.delta_cardinal_cardinal_CI;

deltaI_cardinal_cross       = results_fisher_beforelearning.delta_cardinal_oblique_median;
deltaI_cardinal_cross_CI    = results_fisher_beforelearning.delta_cardinal_oblique_CI;

deltaI_oblique              = results_fisher_beforelearning.delta_oblique_oblique_median;
deltaI_oblique_CI           = results_fisher_beforelearning.delta_oblique_oblique_CI;

deltaI_oblique_cross        = results_fisher_beforelearning.delta_oblique_cardinal_median;
deltaI_oblique_cross_CI     = results_fisher_beforelearning.delta_oblique_cardinal_CI;

fig_it.plot_bar_errorbar(1, deltaI_cardinal, deltaI_cardinal_CI, color_C, color_C)
fig_it.plot_bar_errorbar(3, deltaI_cardinal_cross, deltaI_cardinal_cross_CI, color_C, color_O)
fig_it.plot_bar_errorbar(6, deltaI_oblique, deltaI_oblique_CI, color_O, color_O)
fig_it.plot_bar_errorbar(8, deltaI_oblique_cross, deltaI_oblique_cross_CI, color_O, color_C)

ylim([0, 0.4])
text(0.5,0.38,'Cardinal','color','red','FontSize',14)
text(5.5,0.38,'Oblique','color','blue','FontSize',14)
set(gca,'xtick',[]); set(gca,'fontsize',14)
set(gca,'xtick',[1,3,6,8],'xticklabels',{'Within';'Cross';'Within';'Cross'})
ylabel('$I_\textrm{redundancy}$','Interpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
title('Before learning')

ax_2 = subplot(1,2,2); hold on
set(ax_2,'position',get(ax_2,'position')+[0.01 0.03 0.04 -0.03]);

deltaI_cardinal             = results_fisher_afterlearning.delta_cardinal_cardinal_median;
deltaI_cardinal_CI          = results_fisher_afterlearning.delta_cardinal_cardinal_CI;

deltaI_cardinal_cross       = results_fisher_afterlearning.delta_cardinal_oblique_median;
deltaI_cardinal_cross_CI    = results_fisher_afterlearning.delta_cardinal_oblique_CI;

deltaI_oblique              = results_fisher_afterlearning.delta_oblique_oblique_median;
deltaI_oblique_CI           = results_fisher_afterlearning.delta_oblique_oblique_CI;

deltaI_oblique_cross        = results_fisher_afterlearning.delta_oblique_cardinal_median;
deltaI_oblique_cross_CI     = results_fisher_afterlearning.delta_oblique_cardinal_CI;

fig_it.plot_bar_errorbar(1, deltaI_cardinal, deltaI_cardinal_CI, color_C, color_C)
fig_it.plot_bar_errorbar(3, deltaI_cardinal_cross, deltaI_cardinal_cross_CI, color_C, color_O)
fig_it.plot_bar_errorbar(6, deltaI_oblique, deltaI_oblique_CI, color_O, color_O)
fig_it.plot_bar_errorbar(8, deltaI_oblique_cross, deltaI_oblique_cross_CI, color_O, color_C)

ylim([0, 0.4])
text(0.5,0.38,'Cardinal','color','red','FontSize',14)
text(5.5,0.38,'Oblique','color','blue','FontSize',14)
set(gca,'xtick',[]); set(gca,'fontsize',14)
set(gca,'xtick',[1,3,6,8],'xticklabels',{'Within';'Cross';'Within';'Cross'})
%ylabel('$I_\textrm{redundancy}$','Interpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
title('After learning')

save_folder = '../../figures/figures_final/model_fisher';
save_name = fullfile(save_folder,'model_redundancy_real_cross.svg');
print(save_name,'-dsvg','-vector')
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