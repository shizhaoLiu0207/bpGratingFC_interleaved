clear all
clc
close all
global   bpGlobal  ftsize
bpGratingFCGlobal();
ftsize = 16;
figFolder = '../../figures/figures_informal/Bayesian_model_simulation';
%%
dataFolder = '../../results/neural/fisherInfo_direct/fisherInfo_direct_modelInterleaved_versionControl/subset_32_random_1000';

%%% fisher info results of coherence-combined
dataName = fullfile(dataFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_syntheticData');
load(dataName);
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
%%
nSession = numel(results_all);
CI_level = 68;
prc = [(100-CI_level)/2, 100 - (100-CI_level)/2];
[fisher_cardinal_cross, fisher_oblique_cross] = deal(zeros(nSession, 1));
[fisher_cardinal_cross_CI, fisher_oblique_cross_CI] = deal(zeros(nSession, 2));
[b_PF, cardinal_prior, oblique_prior, cardinal_delta, oblique_delta] = deal(zeros(nSession, 1));
for n = 1:numel(results_all)
    fisher_cardinal_cross_sample = results_all(n).combine_fisher_cardinal_oblique_sample - ...
        results_cross_sizeControl(n).combine_fisher_cardinal_cardinal_sample;
    fisher_oblique_cross_sample = results_all(n).combine_fisher_oblique_cardinal_sample - ...
        results_cross_sizeControl(n).combine_fisher_oblique_oblique_sample;


    fisher_cardinal_cross(n) = median(fisher_cardinal_cross_sample);
    fisher_oblique_cross(n) = median(fisher_oblique_cross_sample);
    fisher_cardinal_cross_CI(n,:) = prctile(fisher_cardinal_cross_sample, prc);
    fisher_oblique_cross_CI(n,:)  = prctile(fisher_oblique_cross_sample, prc);

    tokens = regexp(results_all(n).sessionStr, 'Model_bPF_([-]?[\d_]+)_cardinal_delta_([\d_]+)_prior_([\d_]+)_oblique_delta_([\d_]+)_prior_([\d_]+)', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params            = tokens{1}; % Extract matched tokens
    b_PF_str                    = strrep(extracted_params{1}, '_', '.'); % Replace _ with .
    cardinal_delta_str          = strrep(extracted_params{2}, '_', '.');
    cardinal_prior_str          = strrep(extracted_params{3}, '_', '.');
    oblique_delta_str           = strrep(extracted_params{4}, '_', '.');
    oblique_prior_str           = strrep(extracted_params{5}, '_', '.');
    
    b_PF(n)                     = str2double(b_PF_str);
    cardinal_delta(n)           = str2double(cardinal_delta_str);
    cardinal_prior(n)           = str2double(cardinal_prior_str);
    oblique_delta(n)            = str2double(oblique_delta_str);
    oblique_prior(n)            = str2double(oblique_prior_str);       
end
[~,i_sort] = sort(b_PF,'ascend');
b_PF = b_PF(i_sort);
cardinal_delta = cardinal_delta(i_sort);
cardinal_prior = cardinal_prior(i_sort);
oblique_delta = oblique_delta(i_sort);
oblique_prior = oblique_prior(i_sort);
fisher_cardinal_cross = fisher_cardinal_cross(i_sort);
fisher_oblique_cross   = fisher_oblique_cross(i_sort);
fisher_cardinal_cross_CI = fisher_cardinal_cross_CI(i_sort,:);
fisher_oblique_cross_CI   = fisher_oblique_cross_CI(i_sort,:);

b_PF_list = unique(b_PF); 
prior_list = unique(cardinal_prior);
delta_list = unique(cardinal_delta);
%% 1. b_PF = 0, in a uniform population, check effect of imbalanced delta 
figure
set(gcf,'units','normalized','position',[0,0,1,1]);

for i = 1:numel(delta_list)
    for j = 1:numel(prior_list)
        subplot(numel(delta_list),numel(prior_list),(i-1)*numel(prior_list) + j)
        idx = b_PF == 0 & cardinal_prior == prior_list(j) &  oblique_prior == prior_list(j) & cardinal_delta == delta_list(i);
        
        errorbar([1:sum(idx)],fisher_cardinal_cross(idx),fisher_cardinal_cross(idx) - fisher_cardinal_cross_CI(idx,1),...
            fisher_cardinal_cross_CI(idx,2) - fisher_cardinal_cross(idx), 'linewidth',2,'color',bpGlobal.color_list.color_cardinal);
        hold on
        errorbar([1:sum(idx)],fisher_oblique_cross(idx),fisher_oblique_cross(idx) - fisher_oblique_cross_CI(idx,1),...
            fisher_oblique_cross_CI(idx,2) - fisher_oblique_cross(idx), 'linewidth',2,'color',bpGlobal.color_list.color_oblique);

        line([0.5, sum(idx) + 0.5],[0,0],'color','black','linestyle','--','linewidth',1.5)
        xlim([0.5, sum(idx) + 0.5])
        ylabel('I_{cross} - I_{real}');
        set(gca,'xtick',[1:sum(idx)],'xticklabels',oblique_delta(idx))

        xlabel('Oblique delta')
        title(sprintf('Prior = %.2f, cardinal delta = %.2f', prior_list(j), delta_list(i)))
        set(gca,'fontsize',18)
    
    end
end
sgtitle('Effect of delta, uniform phi_x','fontsize',20)
savename = fullfile(figFolder,'fisherInfo_cross_delta.png');
print(gcf,savename,'-dpng');
%% 2. b_PF = 0, in a uniform population, check effect of prior
figure
set(gcf,'units','normalized','position',[0,0,1,1]);

for i = 1:numel(prior_list)
    for j = 1:numel(delta_list)
       
        subplot(numel(prior_list),numel(delta_list),(i-1)*numel(delta_list) + j)
        idx = b_PF == 0 & cardinal_delta == delta_list(j) &  oblique_delta == delta_list(j) & cardinal_prior == prior_list(i);
        
        errorbar([1:sum(idx)],fisher_cardinal_cross(idx),fisher_cardinal_cross(idx) - fisher_cardinal_cross_CI(idx,1),...
            fisher_cardinal_cross_CI(idx,2) - fisher_cardinal_cross(idx), 'linewidth',2,'color',bpGlobal.color_list.color_cardinal);
        hold on
        errorbar([1:sum(idx)],fisher_oblique_cross(idx),fisher_oblique_cross(idx) - fisher_oblique_cross_CI(idx,1),...
            fisher_oblique_cross_CI(idx,2) - fisher_oblique_cross(idx), 'linewidth',2,'color',bpGlobal.color_list.color_oblique);

        line([0.5, sum(idx) + 0.5],[0,0],'color','black','linestyle','--','linewidth',1.5)
        xlim([0.5, sum(idx) + 0.5])
        ylabel('I_{cross} - I_{real}');
        set(gca,'xtick',[1:sum(idx)],'xticklabels',oblique_prior(idx))

        xlabel('Oblique prior')
        title(sprintf('Delta = %.2f, cardinal prior = %.2f', delta_list(j), prior_list(i)))
        set(gca,'fontsize',18)
    
    end
end
sgtitle('Effect of prior, uniform phi_x','fontsize',20)
savename = fullfile(figFolder,'fisherInfo_cross_prior.png');
print(gcf,savename,'-dpng');
%% 3. check effect of b_PF
figure
set(gcf,'units','normalized','position',[0,0,1,1]);

for i = 1:numel(prior_list)
    for j = 1:numel(delta_list)
        subplot(numel(prior_list),numel(delta_list),(i-1)*numel(delta_list) + j)
        idx = cardinal_prior == prior_list(i) & oblique_prior == prior_list(i) & ...
            cardinal_delta == delta_list(j) & oblique_delta == delta_list(j);
                
        errorbar([1:sum(idx)],fisher_cardinal_cross(idx),fisher_cardinal_cross(idx) - fisher_cardinal_cross_CI(idx,1),...
            fisher_cardinal_cross_CI(idx,2) - fisher_cardinal_cross(idx), 'linewidth',2,'color',bpGlobal.color_list.color_cardinal);
        hold on
        errorbar([1:sum(idx)],fisher_oblique_cross(idx),fisher_oblique_cross(idx) - fisher_oblique_cross_CI(idx,1),...
            fisher_oblique_cross_CI(idx,2) - fisher_oblique_cross(idx), 'linewidth',2,'color',bpGlobal.color_list.color_oblique);

        line([0.5, sum(idx) + 0.5],[0,0],'color','black','linestyle','--','linewidth',1.5)
        xlim([0.5, sum(idx) + 0.5])
        box off
        set(gca,'xtick',[1:sum(idx)],'xticklabels',b_PF(idx))

        xlabel('b PF')
        ylabel('I_{cross} - I_{real}');
        title(sprintf('Delta = %.2f, prior = %.2f', delta_list(j), prior_list(i)))
        set(gca,'fontsize',18)
    end
end
savename = fullfile(figFolder,'fisherInfo_cross_bPF.png');
print(gcf,savename,'-dpng');
%% 4. find sessions whose cardinal cross fisher is significantly larger than zero but oblique close to zero
idx_base = fisher_cardinal_cross_CI(:,1) > 0 & ...
        ((fisher_oblique_cross > 0 & fisher_oblique_cross_CI(:,1) < 0) | ...
        (fisher_oblique_cross < 0 & fisher_oblique_cross_CI(:,2) > 0));
population = 'balanced';
switch population
    case 'cardinal_biased'
        idx = find(idx_base & b_PF > 0);
    case 'oblique_biased'
        idx = find(idx_base & b_PF < 0);
    case 'balanced'
        idx = find(idx_base & b_PF == 0);
end

figure;
set(gcf,'units','normalized','position',[0,0,1,1]);
nCol = ceil(sqrt(numel(idx)));
nRow = ceil(numel(idx) / nCol);
for n = 1:numel(idx)
    subplot(nCol,nRow,n)
     errorbar(1,fisher_cardinal_cross(idx(n)),fisher_cardinal_cross(idx(n)) - fisher_cardinal_cross_CI(idx(n),1),...
            fisher_cardinal_cross_CI(idx(n),2) - fisher_cardinal_cross(idx(n)), 'linewidth',2,'color',bpGlobal.color_list.color_cardinal);
    hold on
    errorbar(1,fisher_oblique_cross(idx(n)),fisher_oblique_cross(idx(n)) - fisher_oblique_cross_CI(idx(n),1),...
        fisher_oblique_cross_CI(idx(n),2) - fisher_oblique_cross(idx(n)), 'linewidth',2,'color',bpGlobal.color_list.color_oblique);

    line([0.5, 1.5],[0,0],'color','black','linestyle','--','linewidth',1.5)
    xlim([0.5, 1.5])
    ylabel('I_{cross} - I_{real}');
    set(gca,'fontsize',14);
    box off
    title({sprintf('b PF = %.1f',b_PF(idx(n)));sprintf('cardinal prior = %.2f, delta = %.2f',cardinal_prior(idx(n)), cardinal_delta(idx(n)));...
        sprintf('oblique prior = %.2f, delta = %.2f',oblique_prior(idx(n)), oblique_delta(idx(n)))})
end
sgtitle(population,'fontsize',20,'interpreter','none')
savename = fullfile(figFolder,sprintf('parameter_setup_asymmetry_%s.png',population));

print(gcf,savename,'-dpng');