clear all
clc
close all
addpath('/Users/liushizhao/projects_local/sampling_decision');
%%%% This scripts test dual-task-switching mode of the prob Bayesian model
%%%% 
%% Parameter set-up
image_task          = 'cardinal';

stimulus_contrast_list = {[15,0],[12,0],[9,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,9],[0,12],[0,15]};
%%%%% use these priors for dual-task swithcing mode
prior_task_dual_list{1} = [1,0];
prior_task_dual_list{2} = [0.8,0.2];
prior_task_dual_list{3} = [0.5,0.5];
prior_task_dual_list{4} = [0.8,0.5];
prior_task_dual_list{5} = [1,1];
prior_task_dual_list{6} = [0,0];
prior_task_dual_list{7} = [0.9,0.9];
prior_task_dual_list{8} = [0.1,0.1];
%%%%% use these priors for single-task switching mode
prior_task_single_list{1} = [1,0];
prior_task_single_list{2} = [0.8,0.2];
prior_task_single_list{3} = [0.5,0.5];

%%%%% Note: We want to compare single and dual task switching mode. 
%%%%% 1. The first scenario, where one task prior is 1, shoud be equivalent
%%%%% 2. Not sure about the second and the third scenario
%%%%% 3. Scenario 4,5,6 are only possible in the dual-mode. Compare 2 and 4
%%%%% would be interesting. Scenario 5 and 6 are extreme cases, let see
%%%%% what happens
savepath = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/test_dual_task';
%% Generate data
doGenerate = 0;
if doGenerate
    
    %%%% Generate dual task data
    nDual = numel(prior_task_dual_list);
    for n = 1:nDual
        
        prior_task_cardinal = prior_task_dual_list{n}(1);
        prior_task_oblique  = prior_task_dual_list{n}(2);
    
        prior_cardinal_str  = strrep(sprintf('%.1f',prior_task_cardinal), '.', '_');
        prior_oblique_str   = strrep(sprintf('%.1f',prior_task_oblique), '.', '_');
    
        savename = fullfile(savepath,sprintf('Dual_image_%s_prior_cardinal_%s_prior_oblique_%s',image_task, prior_cardinal_str, prior_oblique_str));
        for i = 1:numel(stimulus_contrast_list)
            stimulus_contrast = stimulus_contrast_list{i};
            P   = S_Exp_Para('test-dualTask', 'I.stimulus_contrast',stimulus_contrast,'I.image_task',image_task,...
                            'G.prior_task_cardinal', prior_task_cardinal, 'G.prior_task_oblique', prior_task_oblique);
           
            dat_out(i).mode                 = 'dual';
            dat_out(i).prior_task_cardinal  = prior_task_cardinal;
            dat_out(i).prior_task_oblique   = prior_task_oblique;
            dat_out(i).stimulus_contrast    = stimulus_contrast;
            dat_out(i).dat                  = S_Experiment(P);
        end
        save(savename, 'dat_out')
        clear dat_out
    end
    
    %%% Generate single task data
    nSingle = numel(prior_task_single_list);
    for n = 1:nSingle

        prior_task_cardinal = prior_task_single_list{n}(1);
        prior_task_oblique  = prior_task_single_list{n}(2);

        prior_cardinal_str  = strrep(sprintf('%.1f',prior_task_cardinal), '.', '_');
        prior_oblique_str   = strrep(sprintf('%.1f',prior_task_oblique), '.', '_');

        savename = fullfile(savepath,sprintf('Single_image_%s_prior_cardinal_%s_prior_oblique_%s',image_task,prior_cardinal_str, prior_oblique_str));
        for i = 1:numel(stimulus_contrast_list)
            stimulus_contrast = stimulus_contrast_list{i};
            P   = S_Exp_Para('test-interleaved', 'I.stimulus_contrast',stimulus_contrast,'I.image_task',image_task,...
                            'G.prior_task', [prior_task_cardinal, prior_task_oblique]);

            dat_out(i).mode                 = 'dual';
            dat_out(i).prior_task_cardinal  = prior_task_cardinal;
            dat_out(i).prior_task_oblique   = prior_task_oblique;
            dat_out(i).stimulus_contrast    = stimulus_contrast;
            dat_out(i).dat                  = S_Experiment(P);
        end
        save(savename, 'dat_out')
        clear dat_out
    end
end
%% 1. psychometric curve of all conditions
figure
for n = 1:8
    subplot(3,3,n)
    %%%% dual task
    mode = 'dual';
    prior_task = prior_task_dual_list{n};
    
    image_task = 'cardinal';
    [prob_choice_2_cardinal,sem_choice_2_cardinal, contrast_cardinal]  = extract_choice(image_task, prior_task, mode);
    
    image_task = 'oblique';
    [prob_choice_2_oblique,sem_choice_2_oblique,contrast_oblique]    = extract_choice(image_task, prior_task, mode);
    
    h(1) = errorbar(contrast_cardinal, prob_choice_2_cardinal, sem_choice_2_cardinal,'LineWidth',2,'color','red'); hold on
    h(2) = errorbar(contrast_oblique,  prob_choice_2_oblique, sem_choice_2_oblique,'LineWidth',2,'color','blue'); hold on

    %%%% single task
    if n <= 3
         mode = 'single';
        prior_task = prior_task_single_list{n};
        
        image_task = 'cardinal';
        [prob_choice_2_cardinal,sem_choice_2_cardinal, contrast_cardinal]  = extract_choice(image_task, prior_task, mode);
        
        image_task = 'oblique';
        [prob_choice_2_oblique,sem_choice_2_oblique,contrast_oblique]    = extract_choice(image_task, prior_task, mode);
        
        h(1) = errorbar(contrast_cardinal, prob_choice_2_cardinal, sem_choice_2_cardinal,'LineWidth',2,'color','red','LineStyle','--'); hold on
        h(2) = errorbar(contrast_oblique,  prob_choice_2_oblique, sem_choice_2_oblique,'LineWidth',2,'color','blue','LineStyle','--'); hold on

    end
    set(gca,'fontsize',14)
    xlabel('Contrast')
    ylabel('Prob. Ori 2')
    title(sprintf('Prior_{cardinal} = %.1f, Prior_{oblique} = %.1f',prior_task(1),prior_task(2)));
    grid on
    grid minor
end
%% 2. check variables 
doThis = 0;
if doThis
    use_list = [-12, -6, 0, 6, 12]; % five contrast levels
    savepath = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/test_dual_task';
    %%%% 2.1 T_cardinal = 1, T_oblique = 0 for both mode
    f = figure;
    set(gcf,'unit','normalized','Position',[0,0,1,1]);
    
    single_name = fullfile(savepath, 'Single_image_cardinal_prior_cardinal_1_0_prior_oblique_0_0');
    load(single_name)
    check_var_behavior(dat_out,use_list,[0,0.47,0.5,0.5],f)
    
    cond2_cardinal_name = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_1_0_prior_oblique_0_0');
    load(cond2_cardinal_name)
    check_var_behavior(dat_out,use_list,[0.5,0.47,0.5,0.5],f)
    
    single_name = fullfile(savepath, 'Single_image_oblique_prior_cardinal_1_0_prior_oblique_0_0');
    load(single_name)
    check_var_behavior(dat_out,use_list,[0,0,0.5,0.5],f)
    
    cond2_cardinal_name = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_1_0_prior_oblique_0_0');
    load(cond2_cardinal_name)
    check_var_behavior(dat_out,use_list,[0.5,0,0.5,0.5],f)
    
    delete(findall(gcf,'type','annotation'))
    annotation('textbox',[0.01,0.93,0.15,0.04],'string','Cardinal','fontsize',20,'FontWeight','bold','EdgeColor','none','Color','red')
    annotation('textbox',[0.01,0.46,0.15,0.04],'string','Oblique','fontsize',20,'FontWeight','bold','EdgeColor','none','Color','blue')
    annotation('textbox',[0.2,0.93,0.20,0.04],'string','Single, T_{cardinal} = 1, T_{oblique} = 0','fontsize',16,'FontWeight','bold','EdgeColor','none')
    annotation('textbox',[0.7,0.93,0.20,0.04],'string','Dual, T_{cardinal} = 1, T_{oblique} = 0','fontsize',16,'FontWeight','bold','EdgeColor','none')
    %%%%% 2.2 T_cardinal = 0.8, T_oblique = 0.2 for both mode
    f = figure;
    set(gcf,'unit','normalized','Position',[0,0,1,1]);
    
    single_name = fullfile(savepath, 'Single_image_cardinal_prior_cardinal_0_8_prior_oblique_0_2');
    load(single_name)
    check_var_behavior(dat_out,use_list,[0,0.47,0.5,0.5],f)
    
    cond2_cardinal_name = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_0_8_prior_oblique_0_2');
    load(cond2_cardinal_name)
    check_var_behavior(dat_out,use_list,[0.5,0.47,0.5,0.5],f)
    
    single_name = fullfile(savepath, 'Single_image_oblique_prior_cardinal_0_8_prior_oblique_0_2');
    load(single_name)
    check_var_behavior(dat_out,use_list,[0,0,0.5,0.5],f)
    
    cond2_cardinal_name = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_0_8_prior_oblique_0_2');
    load(cond2_cardinal_name)
    check_var_behavior(dat_out,use_list,[0.5,0,0.5,0.5],f)
    
    delete(findall(gcf,'type','annotation'))
    annotation('textbox',[0.01,0.93,0.15,0.04],'string','Cardinal','fontsize',20,'FontWeight','bold','EdgeColor','none','Color','red')
    annotation('textbox',[0.01,0.46,0.15,0.04],'string','Oblique','fontsize',20,'FontWeight','bold','EdgeColor','none','Color','blue')
    annotation('textbox',[0.2,0.93,0.20,0.04],'string','Single, T_{cardinal} = 0.8, T_{oblique} = 0.2','fontsize',16,'FontWeight','bold','EdgeColor','none')
    annotation('textbox',[0.7,0.93,0.20,0.04],'string','Dual, T_{cardinal} = 0.8, T_{oblique} = 0.2','fontsize',16,'FontWeight','bold','EdgeColor','none')
    
    %%%% 2.3 T_cardinal = 0.8, T_oblique = 0.2/T_oblique = 0.5 for the dual mode
    f = figure;
    set(gcf,'unit','normalized','Position',[0,0,1,1]);
    
    single_name = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_0_8_prior_oblique_0_2');
    load(single_name)
    check_var_behavior(dat_out,use_list,[0,0.47,0.5,0.5],f)
    
    cond2_cardinal_name = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_0_8_prior_oblique_0_5');
    load(cond2_cardinal_name)
    check_var_behavior(dat_out,use_list,[0.5,0.47,0.5,0.5],f)
    
    single_name = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_0_8_prior_oblique_0_2');
    load(single_name)
    check_var_behavior(dat_out,use_list,[0,0,0.5,0.5],f)
    
    cond2_cardinal_name = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_0_8_prior_oblique_0_5');
    load(cond2_cardinal_name)
    check_var_behavior(dat_out,use_list,[0.5,0,0.5,0.5],f)
    
    delete(findall(gcf,'type','annotation'))
    annotation('textbox',[0.01,0.93,0.15,0.04],'string','Cardinal','fontsize',20,'FontWeight','bold','EdgeColor','none','Color','red')
    annotation('textbox',[0.01,0.46,0.15,0.04],'string','Oblique','fontsize',20,'FontWeight','bold','EdgeColor','none','Color','blue')
    annotation('textbox',[0.2,0.93,0.20,0.04],'string','Single, T_{cardinal} = 0.8, T_{oblique} = 0.2','fontsize',16,'FontWeight','bold','EdgeColor','none')
    annotation('textbox',[0.7,0.93,0.20,0.04],'string','Dual, T_{cardinal} = 0.8, T_{oblique} = 0.5','fontsize',16,'FontWeight','bold','EdgeColor','none')
end
%% 3. Directly compare a few variables between conditions
figure
compare_case = 'twoDual';
switch compare_case
    case 'twoDual'
        cond1_cardinal_name     = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_0_8_prior_oblique_0_2');
        cond2_cardinal_name     = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_0_8_prior_oblique_0_5');

        cond1_oblique_name     = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_0_8_prior_oblique_0_2');
        cond2_oblique_name     = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_0_8_prior_oblique_0_5');

        cond_1_str = 'prior_{oblique} = 0.2';
        cond_2_str = 'prior_{oblique} = 0.5';
        sgtitle_str = 'Two dual conditions, prior_{cardinal} = 0.8';

    case 'twoIdeal'
        cond1_cardinal_name       = fullfile(savepath, 'Single_image_cardinal_prior_cardinal_1_0_prior_oblique_0_0.mat');
        cond2_cardinal_name       = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_1_0_prior_oblique_0_0.mat');

        cond1_oblique_name       = fullfile(savepath, 'Single_image_oblique_prior_cardinal_1_0_prior_oblique_0_0.mat');
        cond2_oblique_name       = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_1_0_prior_oblique_0_0.mat');

        cond_1_str = 'Single';
        cond_2_str = 'Dual';
        sgtitle_str = 'prior_{cardinal} = 1.0, prior_{oblique} = 0';
    case 'twoMix'
        cond1_cardinal_name       = fullfile(savepath, 'Single_image_cardinal_prior_cardinal_0_8_prior_oblique_0_2.mat');
        cond2_cardinal_name       = fullfile(savepath, 'Dual_image_cardinal_prior_cardinal_0_8_prior_oblique_0_2.mat');

        cond1_oblique_name       = fullfile(savepath, 'Single_image_oblique_prior_cardinal_0_8_prior_oblique_0_2.mat');
        cond2_oblique_name       = fullfile(savepath, 'Dual_image_oblique_prior_cardinal_0_8_prior_oblique_0_2.mat');

        cond_1_str = 'Single';
        cond_2_str = 'Dual';
        sgtitle_str = 'prior_{cardinal} = 0.8, prior_{oblique} = 0.2';
end

dat_cond1   = load(cond1_cardinal_name);
dat_cond2   = load(cond2_cardinal_name);


i = 4; % low contrast
subplot(3,2,1)
plot(squeeze(mean(dat_cond1.dat_out(i).dat.O(:,2,:),1))); hold on
plot(squeeze(mean(dat_cond2.dat_out(i).dat.O(:,2,:),1)));
ylim([0.4, 1])
ylabel('P(O1) low contrast')
legend(cond_1_str, cond_2_str)
title('Cardinal task')
set(gca,'fontsize',16)

subplot(3,2,3)
phi_g = dat_cond1.dat_out(1).dat.Projection.phi_g;
plot(phi_g, squeeze(mean(mean(dat_cond1.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
plot(phi_g, squeeze(mean(mean(dat_cond2.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
ylabel('P(G) low contrast')
set(gca,'fontsize',16)

i = 6; % zero 
subplot(3,2,5)
phi_g = dat_cond1.dat_out(1).dat.Projection.phi_g;
plot(phi_g, squeeze(mean(mean(dat_cond1.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
plot(phi_g, squeeze(mean(mean(dat_cond2.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
ylabel('P(G) zero contrast')
set(gca,'fontsize',16)


dat_cond1   = load(cond1_oblique_name);
dat_cond2   = load(cond2_oblique_name);



i = 4; 
subplot(3,2,2)
plot(squeeze(mean(dat_cond1.dat_out(i).dat.O(:,2,:),1))); hold on
plot(squeeze(mean(dat_cond2.dat_out(i).dat.O(:,2,:),1)));
ylim([0.4, 1])
ylabel('P(O1) low contrast')
title('Oblique task')
set(gca,'fontsize',16)

subplot(3,2,4)
phi_g = dat_cond1.dat_out(1).dat.Projection.phi_g;
plot(phi_g, squeeze(mean(mean(dat_cond1.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
plot(phi_g, squeeze(mean(mean(dat_cond2.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
ylabel('P(G) low contrast')
set(gca,'fontsize',16)

i = 6;
subplot(3,2,6)
phi_g = dat_cond1.dat_out(1).dat.Projection.phi_g;
plot(phi_g, squeeze(mean(mean(dat_cond1.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
plot(phi_g, squeeze(mean(mean(dat_cond2.dat_out(i).dat.G(:,:,:,20:100),4),1))); hold on
ylabel('P(G) zero contrast')
set(gca,'fontsize',16)

sgtitle(sgtitle_str,'fontsize',18);
%% help functions

function check_var_behavior(dat_out,use_list, pos, parent)
for i = 1:11
    stimulus_contrast_list(i) = dat_out(i).stimulus_contrast(2) - dat_out(i).stimulus_contrast(1);
end
color_list = {[0, 0.4470, 0.7410], [0.9290, 0.6940, 0.1250] ,[0.25, 0.25, 0.25],...
            [0.9290, 0.6940, 0.1250],[0, 0.4470, 0.7410]};
style_list = {'-','-',':','--','--';};

if nargin < 1 || isempty(pos)
    pos = [0 0 1 1];   % full figure
end
if nargin < 2 || isempty(parent)
    parent = gcf;
end

% layout: 3 rows x 2 columns
nRows = 3;
nCols = 2;

% margins (in normalized units within pos)
left   = 0.06;
right  = 0.02;
bottom = 0.08;
top    = 0.04;
hgap   = 0.04;   % horiz gap between subplots
vgap   = 0.06;   % vert gap

innerW = pos(3) - left - right;
innerH = pos(4) - top  - bottom;

tileW  = (innerW - (nCols-1)*hgap) / nCols;
tileH  = (innerH - (nRows-1)*vgap) / nRows;

ax = gobjects(nRows*nCols,1);
k = 0;

for r = 1:nRows
    for c = 1:nCols
        k = k + 1;

        x = pos(1) + left + (c-1)*(tileW + hgap);
        % top row has largest y; convert row index to y position
        y = pos(2) + bottom + (nRows - r)*(tileH + vgap);

        ax(k) = axes('Parent', parent, ...
                     'Position', [x y tileW tileH]);
        
        if r == 1 & c == 1
            for i = 1:numel(use_list)
                idx = stimulus_contrast_list == use_list(i);
                if isfield(dat_out(idx).dat, 'T_cardinal')
                    plot(squeeze(mean(dat_out(idx).dat.T_cardinal(:,2,:),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2); % p_T_cardinal
                else
                    plot(squeeze(mean(dat_out(idx).dat.T(:,2,:),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2);
                end
                hold on
            end
            ylim([-0.05,1.05]); ylabel('P(T_{cardinal})');xlabel('Sample'); set(gca,'fontsize',14);
        end

        if r == 1 & c == 2
            for i = 1:numel(use_list)
                idx = stimulus_contrast_list == use_list(i);
                if isfield(dat_out(idx).dat, 'T_oblique')
                    plot(squeeze(mean(dat_out(idx).dat.T_oblique(:,2,:),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2); % p_T_oblique
                else
                    plot(squeeze(mean(dat_out(idx).dat.T(:,3,:),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2);
                end
                
                hold on
            
            end
            ylim([-0.05,1.05]);ylabel('P(T_{oblique})');xlabel('Sample'); set(gca,'fontsize',14);
        end

        if r == 2 & c == 1 
            for i = 1:numel(use_list)
                idx = stimulus_contrast_list == use_list(i);
                plot(squeeze(mean(dat_out(idx).dat.O(:,2,:),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2); hold on
               % plot(squeeze(mean(dat_out(i).dat.O(:,3,:),1))); 
            end
            ylim([-0.05,1.05]);ylabel('P(O1)');xlabel('Sample'); set(gca,'fontsize',14);
        end

        if r == 2 & c == 2  
            for i = 1:numel(use_list)
                idx = stimulus_contrast_list == use_list(i);
                plot(squeeze(mean(dat_out(idx).dat.O(:,3,:),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2); hold on
               % plot(squeeze(mean(dat_out(i).dat.O(:,3,:),1))); 
            end
            ylim([-0.05,1.05]);ylabel('P(O2)');xlabel('Sample'); set(gca,'fontsize',14);
        end
        
        if r == 3 & c == 1 
            phi_g = dat_out(1).dat.Projection.phi_g;
            
            for i = 1:numel(use_list)
                idx = stimulus_contrast_list == use_list(i);
                plot(phi_g, squeeze(mean(mean(dat_out(idx).dat.G(:,:,:,20:40),4),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2); hold on
            end
            ylim([-0.05,0.65]); ylabel('P(G)');xlabel('phi_g'); set(gca,'fontsize',14);
        end

        if r == 3 & c == 2  
            phi_g = dat_out(1).dat.Projection.phi_g;
            for i = 1:numel(use_list)
                idx = stimulus_contrast_list == use_list(i);
                plot(phi_g, squeeze(mean(mean(dat_out(idx).dat.G(:,:,:,80:100),4),1)),'Color',color_list{i},'linestyle',style_list{i},'linewidth',2); hold on
            end
            ylim([-0.05,0.65]); ylabel('P(G)');xlabel('phi_g'); set(gca,'fontsize',14);
        end
    end
end
end
function [prob_choice_2,sem_choice_2, contrast] = extract_choice(image_task, prior_task, mode)
    savepath = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/test_dual_task';
    prior_task_cardinal = prior_task(1);
    prior_task_oblique  = prior_task(2);
    prior_cardinal_str  = strrep(sprintf('%.1f',prior_task_cardinal), '.', '_');
    prior_oblique_str   = strrep(sprintf('%.1f',prior_task_oblique), '.', '_');
    
    savename = fullfile(savepath,sprintf('%s_image_%s_prior_cardinal_%s_prior_oblique_%s',mode, image_task,prior_cardinal_str, prior_oblique_str));
    D = load(savename);

    N = numel(D.dat_out);
    [contrast,sem_choice_2, prob_choice_2] = deal(zeros(N,1));

    for n = 1:N
        contrast(n)         = D.dat_out(n).stimulus_contrast(2) - D.dat_out(n).stimulus_contrast(1);
        decision            = (D.dat_out(n).dat.O(:,3,end) > 0.5) + 1;
        prob_choice_2(n)    = sum(decision == 2) / numel(decision);
        sem_choice_2(n)     = sqrt( prob_choice_2(n) * (1 - prob_choice_2(n)) / numel(decision));
    end
    
end

