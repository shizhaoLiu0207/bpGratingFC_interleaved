clear all
clc
close all
global bpGlobal
bpGratingFCGlobal();

%%%% This script aims to check psychometric kernel of the model with the
%%%% parameters that best explain the neural data.


%%%% monkey R: delta_cardinal = 0.06, prior_cardinal = 0.60, 
%%%% delta_oblique = 0.04, prior_oblique = 0.80;
%%%% monkey G: delta_cardinal = 0.08, prior_cardinal = 1.00,
%%%% delta_oblique = 0.04, prior_oblique = 0.60
monkeyR_delta_cardinal = 0.06;
monkeyR_prior_cardinal = 0.60;
monkeyR_delta_oblique = 0.04;
monkeyR_prior_oblique = 0.80;

monkeyG_delta_cardinal = 0.08;
monkeyG_prior_cardinal = 1.00;
monkeyG_delta_oblique = 0.04;
monkeyG_prior_oblique = 0.60;

%% Generate synthetic data with pre-set parameters

doGenerate = 0;
%%% on macbook
home_folder                     = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel';
%%%% on linux
if ~exist(home_folder)
    home_folder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/data_for_psykernel';
end


if doGenerate




runOptions.number_samples_per_evidence         = 6;
runOptions.stimulus_regime                     = 'dynamic-switching-signal-blocked';

runOptions.delta_list                          = [monkeyR_delta_cardinal, monkeyR_delta_oblique, ...
                                        monkeyG_delta_cardinal, monkeyG_delta_oblique];
runOptions.prior_task_list                     = {[monkeyR_prior_cardinal, 1 - monkeyR_prior_cardinal];...
                                       [1 - monkeyR_prior_oblique, monkeyR_prior_oblique];...
                                       [monkeyG_prior_cardinal, 1 - monkeyG_prior_cardinal];...
                                       [1 - monkeyG_prior_oblique, monkeyG_prior_oblique]};

runOptions.image_task_list                     = {'cardinal';'oblique';'cardinal';'oblique'};


runOptions.nRepeats                            = 5000;

runOptions.stimulus_contrast_list              = {[0,6],[6,0]};
runOptions.nRepeats_list                       = ones(size(runOptions.stimulus_contrast_list)) * runOptions.nRepeats;

runOptions.run_ori_energy           = true;
runOptions.n_ori_bin                = 12; 
runOptions.task_mode                = 'single';
%%%%%%%% 01/06/2025, change back to non-clamp prior
runOptions.clamp_prior              = false;
runOptions.save_folder              = fullfile(home_folder, 'best_parameters_contrast6');
runOptions.nSession                 = 5; %%% run 5 sessions per condition

generate_data_model_psyKernel(runOptions);

end
%% Compute psychometric kernel
doCompute = 0;
version_name = 'best_parameters_contrast6';
data_folder = fullfile(home_folder, version_name);
energy_use = 'proj';
save_name = sprintf('../../results/model_behavior/psyKernel_model_%s_energy_%s',version_name, energy_use);

if doCompute
    
    data_list = dir(fullfile(data_folder,'*.mat'));
    data_list(contains({data_list(:).name},'raw_data')) = [];
    
    psyKernel_model = struct();
    for n = 1:numel(data_list)
        fprintf('Computing session %d/%d \n',n,numel(data_list));
        load(fullfile(data_folder, data_list(n).name));
        
        [map_params,orientationsDEG] = util_it.compute_psyKernel_model(data_use, energy_use);
      
        tokens = regexp(data_list(n).name, '([a-zA-Z0-9]+)_image_task_([a-zA-Z0-9]+)_delta_([-]?[\d_]+)_prior_([\d_]+)_session_([\d_]+)', 'tokens');
        % Convert extracted strings back to numbers
        extracted_params        = tokens{1}; % Extract matched tokens
        task_mode_str           = extracted_params{1};
        imagetask_str           = extracted_params{2};
        delta_str               = strrep(extracted_params{3}, '_', '.');
        prior_str               = strrep(extracted_params{4}, '_', '.');
        session_num_str         = extracted_params{5};

        psyKernel_model(n).mode             = task_mode_str;
        psyKernel_model(n).image_task       = imagetask_str;
        psyKernel_model(n).delta            = str2double(delta_str);
        psyKernel_model(n).prior            = str2double(prior_str);
        psyKernel_model(n).sessionNum       = str2double(session_num_str);
        psyKernel_model(n).w_ori            = map_params.w_ori;
        psyKernel_model(n).w_time           = map_params.w_time;
        psyKernel_model(n).orientationsDEG  = orientationsDEG;
    
        psyKernel_model(n).energy_use       = energy_use;

        %%%% number of repeats and list of coherence
        psyKernel_model(n).cohr_list        = unique(abs(cell2mat({data_use(:).contrast})));
        for k = 1:numel(psyKernel_model(n).cohr_list)
            idx = cell2mat({data_use(:).contrast}) == psyKernel_model(n).cohr_list;
            psyKernel_model(n).num_repeats(k) = numel(data_use(idx).decision);
        end
    end
    save(save_name,'psyKernel_model');
else
    load(save_name);
end


%% visualization shape of the kernel
subjectCode = 'Gr';
kernel_type             = 'spatial';
switch subjectCode
    case 'Ro'
        prior_cardinal_plot = monkeyR_prior_cardinal;
        delta_cardinal_plot = monkeyR_delta_cardinal;
        prior_oblique_plot = monkeyR_prior_oblique;
        delta_oblique_plot = monkeyR_delta_oblique;
    case 'Gr'
        prior_cardinal_plot = monkeyG_prior_cardinal;
        delta_cardinal_plot = monkeyG_delta_cardinal;
        prior_oblique_plot = monkeyG_prior_oblique;
        delta_oblique_plot = monkeyG_delta_oblique;

end
figure
subplot(1,2,1)
plotOptions.mode            = 'Single';
plotOptions.image_task      = 'cardinal';
plotOptions.prior_task      = prior_cardinal_plot;
plotOptions.delta           = delta_cardinal_plot;
plotOptions.plotIndividual  = true;
[r_cardinal, r_oblique,~, w_ori_cardinal_model] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, plotOptions);
title(sprintf('Delta = %.2f,  Prior = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    plotOptions.delta, plotOptions.prior_task, r_cardinal, r_oblique));

subplot(1,2,2)
plotOptions.mode            = 'Single';
plotOptions.image_task      = 'oblique';
plotOptions.prior_task      = prior_oblique_plot;
plotOptions.delta           = delta_oblique_plot;
plotOptions.plotIndividual  = true;
[r_cardinal, r_oblique,~, w_ori_oblique_model] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, plotOptions);
title(sprintf('Delta = %.2f,  Prior = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    plotOptions.delta, plotOptions.prior_task, r_cardinal, r_oblique));

sgtitle(sprintf('Using best parameters for %s',subjectCode),'fontsize',20,'fontweight','bold')
%% directly compare to the empirical kernels
psyKernel_Ro = load('/Users/liushizhao/projects_local/bpGratingEx/results/behavior/Rolo_psyKernel_table_final');
psyKernel_Gr = load('/Users/liushizhao/projects_local/bpGratingEx/results/behavior/Gremlin_psyKernel_table');
switch subjectCode
    case 'Ro'
        sessionName_list  = bpGlobal.rolo.session_list.switching;
        idx_list                        = ismember({psyKernel_Ro.psyKernel_table(:).sessionName},sessionName_list);
        
        tmp = {psyKernel_Ro.psyKernel_table(idx_list).w_ori_cardinal};
        w_ori_cardinal_empirical = cat(2,tmp{:});
        tmp = {psyKernel_Ro.psyKernel_table(idx_list).w_ori_oblique};
        w_ori_oblique_empirical = cat(2,tmp{:});
    case 'Gr'
        sessionName_list  = bpGlobal.gremlin.session_list.interleaved_real;
        idx_list                        = ismember({psyKernel_Gr.psyKernel_table(:).sessionName},sessionName_list);
        
        tmp = {psyKernel_Gr.psyKernel_table(idx_list).w_ori_cardinal};
        w_ori_cardinal_empirical = cat(2,tmp{:});
        tmp = {psyKernel_Gr.psyKernel_table(idx_list).w_ori_oblique};
        w_ori_oblique_empirical = cat(2,tmp{:});
end



nModel = size(w_ori_cardinal_model,2);
nEmpirical = size(w_ori_cardinal_empirical,2);

orientationsDEG = psyKernel_model(1).orientationsDEG;
ideal_kernel_cardinal            = 0.4*sin(pi * (orientationsDEG - 45)/90)';
ideal_kernel_oblique             = 0.4*sin(pi * (orientationsDEG - 90)/90)';

r_model_empirical_cardinal_avg = corr(mean(w_ori_cardinal_model,2),mean(w_ori_cardinal_empirical,2));
r_model_empirical_oblique_avg  = corr(mean(w_ori_oblique_model,2),mean(w_ori_oblique_empirical,2));

r_model_ideal_cardinal_cardinal_avg = corr(mean(w_ori_cardinal_model,2),ideal_kernel_cardinal);
r_model_ideal_cardinal_oblique_avg  = corr(mean(w_ori_cardinal_model,2), ideal_kernel_oblique);

r_model_ideal_oblique_cardinal_avg = corr(mean(w_ori_oblique_model,2),ideal_kernel_cardinal);
r_model_ideal_oblique_oblique_avg  = corr(mean(w_ori_oblique_model,2), ideal_kernel_oblique);


r_empirical_ideal_cardinal_cardinal_avg = corr(mean(w_ori_cardinal_empirical,2),ideal_kernel_cardinal);
r_empirical_ideal_cardinal_oblique_avg  = corr(mean(w_ori_cardinal_empirical,2), ideal_kernel_oblique);

r_empirical_ideal_oblique_cardinal_avg = corr(mean(w_ori_oblique_empirical,2),ideal_kernel_cardinal);
r_empirical_ideal_oblique_oblique_avg  = corr(mean(w_ori_oblique_empirical,2), ideal_kernel_oblique);

r_model_empirical_cardinal_all = zeros(nModel,nEmpirical);
r_model_empirical_oblique_all  = zeros(nModel,nEmpirical);

r_model_ideal_cardinal_cardinal_all = zeros(nModel,1);
r_model_ideal_cardinal_oblique_all  = zeros(nModel,1);
r_model_ideal_oblique_cardinal_all = zeros(nModel,1);
r_model_ideal_oblique_oblique_all  = zeros(nModel,1);

r_empirical_ideal_cardinal_cardinal_all = zeros(nEmpirical,1);
r_empirical_ideal_cardinal_oblique_all  = zeros(nEmpirical,1);
r_empirical_ideal_oblique_cardinal_all = zeros(nEmpirical,1);
r_empirical_ideal_oblique_oblique_all  = zeros(nEmpirical,1);

for i = 1:nModel
    r_model_ideal_cardinal_cardinal_all(i)  = corr(w_ori_cardinal_model(:,i), ideal_kernel_cardinal);
    r_model_ideal_cardinal_oblique_all(i)   = corr(w_ori_cardinal_model(:,i), ideal_kernel_oblique);
    r_model_ideal_oblique_cardinal_all(i)   = corr(w_ori_oblique_model(:,i), ideal_kernel_cardinal);
    r_model_ideal_oblique_oblique_all(i)    = corr(w_ori_oblique_model(:,i), ideal_kernel_oblique);
    for j = 1:nEmpirical
        r_empirical_ideal_cardinal_cardinal_all(j) = corr(w_ori_cardinal_empirical(:,j), ideal_kernel_cardinal);
        r_empirical_ideal_cardinal_oblique_all(j) = corr(w_ori_cardinal_empirical(:,j), ideal_kernel_oblique);
        r_empirical_ideal_oblique_cardinal_all(j) = corr(w_ori_oblique_empirical(:,j), ideal_kernel_cardinal);
        r_empirical_ideal_oblique_oblique_all(j) = corr(w_ori_oblique_empirical(:,j), ideal_kernel_oblique);

        r_model_empirical_cardinal_all(i,j) = corr(w_ori_cardinal_model(:,i), w_ori_cardinal_empirical(:,j));
        r_model_empirical_oblique_all(i,j) = corr(w_ori_oblique_model(:,i), w_ori_oblique_empirical(:,j));
    end
end


edge = [-1:0.1:1];
figure
subplot(3,2,1)
histogram(r_model_empirical_cardinal_all(:),edge, 'facecolor','red','FaceAlpha',0.5);
hold on
line([r_model_empirical_cardinal_avg, r_model_empirical_cardinal_avg], [0,25],'color','black','linewidth',3)
title('corr-model-empirical, cardinal')
set(gca,'fontsize',16)

subplot(3,2,2)
histogram(r_model_empirical_oblique_all(:),edge,'facecolor','blue','FaceAlpha',0.5);
hold on
line([r_model_empirical_oblique_avg, r_model_empirical_oblique_avg], [0,25],'color','black','linewidth',3)
title('corr-model-empirical, oblique')
set(gca,'fontsize',16)

subplot(3,2,3)
histogram(r_model_ideal_cardinal_cardinal_all(:),edge,'facecolor','red','FaceAlpha',0.5);
hold on
h(1) = line([r_model_ideal_cardinal_cardinal_avg, r_model_ideal_cardinal_cardinal_avg], [0,10],...
    'color','black','linewidth',3);
histogram(r_empirical_ideal_cardinal_cardinal_all(:),edge,'facecolor','red','FaceAlpha',0.2);
hold on
h(2) = line([r_empirical_ideal_cardinal_cardinal_avg, r_empirical_ideal_cardinal_cardinal_avg], [0,10],...
    'color','black','linewidth',3,'linestyle','--');
title('corr w. ideal, cardinal')
set(gca,'fontsize',16)
legend(h, 'Model', 'Empirical')


subplot(3,2,4)
histogram(r_model_ideal_oblique_oblique_all(:),edge,'facecolor','blue','FaceAlpha',0.5);
hold on
h(1) = line([r_model_ideal_oblique_oblique_avg, r_model_ideal_oblique_oblique_avg], ...
    [0,10],'color','black','linewidth',3);
histogram(r_empirical_ideal_oblique_oblique_all(:),edge,'facecolor','blue','FaceAlpha',0.2);
hold on
h(2) = line([r_empirical_ideal_oblique_oblique_avg, r_empirical_ideal_oblique_oblique_avg], [0,10],...
    'color','black','linewidth',3,'linestyle','--');
title('corr w. ideal, oblique')
set(gca,'fontsize',16)
legend(h, 'Model', 'Empirical')

subplot(3,2,5)
histogram(r_model_ideal_cardinal_oblique_all(:),edge,'facecolor','red','FaceAlpha',0.5);
hold on
h(1) = line([r_model_ideal_cardinal_oblique_avg, r_model_ideal_cardinal_oblique_avg], [0,10], ...
        'color','black','linewidth',3);
histogram(r_empirical_ideal_cardinal_oblique_all(:),edge,'facecolor','red','FaceAlpha',0.2);
hold on
h(2) = line([r_empirical_ideal_cardinal_oblique_avg, r_empirical_ideal_cardinal_oblique_avg], [0,10], ...
    'color','black','linewidth',3,'linestyle','--');
title('corr w. ideal (cross), cardinal')
set(gca,'fontsize',16)
legend(h, 'Model', 'Empirical')

subplot(3,2,6)
histogram(r_model_ideal_oblique_cardinal_all(:),edge,'facecolor','blue','FaceAlpha',0.5);
hold on
h(1) = line([r_model_ideal_oblique_cardinal_avg, r_model_ideal_oblique_cardinal_avg], [0,10], ...
    'color','black','linewidth',3);
histogram(r_empirical_ideal_oblique_cardinal_all(:),edge,'facecolor','blue','FaceAlpha',0.2);
hold on
h(2) = line([r_empirical_ideal_oblique_cardinal_avg, r_empirical_ideal_oblique_cardinal_avg], [0,10], ...
    'color','black','linewidth',3,'linestyle','--');
title('corr w. ideal (cross), oblique')
set(gca,'fontsize',16)
legend(h, 'Model', 'Empirical')

sgtitle(sprintf('Using best parameters for %s',subjectCode),'fontsize',20,'fontweight','bold')

%% cross-prediction 

doPred = 0;
save_folder_pred = '../../results/model_behavior/predCross/bestparams';
if doPred 
    nFold = 20;
    subjectCode = 'Gr';
    %predResults = run_cross_predict_model(dat_cardinal,dat_oblique, nFold);
    
    
    data_list = dir(fullfile(data_folder,'*.mat'));
    data_list(contains({data_list(:).name},'raw_data')) = [];
    data_list(contains({data_list(:).name},'predCross')) = [];
    data_info = struct();
    for n = 1:numel(data_list)
      
        tokens = regexp(data_list(n).name, '([a-zA-Z0-9]+)_image_task_([a-zA-Z0-9]+)_delta_([-]?[\d_]+)_prior_([\d_]+)_session_([\d_]+)', 'tokens');
        % Convert extracted strings back to numbers
        extracted_params        = tokens{1}; % Extract matched tokens
        task_mode_str           = extracted_params{1};
        imagetask_str           = extracted_params{2};
        delta_str               = strrep(extracted_params{3}, '_', '.');
        prior_str               = strrep(extracted_params{4}, '_', '.');
        session_num_str         = extracted_params{5};
    
    
    
        data_info(n).filename           =  data_list(n).name;
        data_info(n).mode               = task_mode_str;
        data_info(n).image_task         = imagetask_str;
        data_info(n).delta              = str2double(delta_str);
        data_info(n).prior              = str2double(prior_str);
        data_info(n).sessionNum         = str2double(session_num_str);
    end
    
    switch subjectCode
        case 'Ro'
    
            idx_cardinal = find(strcmp({data_info(:).image_task}, 'cardinal') &...
                        [data_info(:).delta] == monkeyR_delta_cardinal &...
                        [data_info(:).prior] == monkeyR_prior_cardinal);
            
            idx_oblique = find(strcmp({data_info(:).image_task}, 'oblique') &...
                        [data_info(:).delta] == monkeyR_delta_oblique &...
                        [data_info(:).prior] == monkeyR_prior_oblique);
        case 'Gr'
            idx_cardinal = find(strcmp({data_info(:).image_task}, 'cardinal') &...
                        [data_info(:).delta] == monkeyG_delta_cardinal &...
                        [data_info(:).prior] == monkeyG_prior_cardinal);
            
            idx_oblique = find(strcmp({data_info(:).image_task}, 'oblique') &...
                        [data_info(:).delta] == monkeyG_delta_oblique &...
                        [data_info(:).prior] == monkeyG_prior_oblique);
    
    end
    
    [A, B] = ndgrid(idx_cardinal, idx_oblique);
    allComb = [A(:) B(:)];
    
    nCombination = size(allComb,1);
    for k = 1:nCombination
        save_name           = fullfile(save_folder_pred,sprintf('predCross_bestparams_monkey_%s_combination_%d.mat',...
                                   subjectCode, k)); 
        if isfile(save_name)
            continue
        end
        fprintf('Running combination %d/%d \n',k, nCombination );
        file_name_cardinal  = fullfile(data_folder, data_info(allComb(k,1)).filename);
        file_name_oblique   = fullfile(data_folder, data_info(allComb(k,2)).filename);
        
        dat_cardinal        = load(file_name_cardinal);
        dat_oblique         = load(file_name_oblique);
        predResults         = run_cross_predict_model(dat_cardinal,dat_oblique, nFold);
        
        predCross_results   = struct();
        predCross_results.mode_cardinal_session          = data_info(allComb(k,1)).mode;
        predCross_results.mode_oblique_session           = data_info(allComb(k,2)).mode;
        predCross_results.delta_cardinal_session         = data_info(allComb(k,1)).delta;
        predCross_results.delta_oblique_session          = data_info(allComb(k,2)).delta;
        predCross_results.prior_cardinal_session         = data_info(allComb(k,1)).prior;
        predCross_results.prior_oblique_session          = data_info(allComb(k,2)).prior;
        predCross_results.sessionNum_cardinal_session    = data_info(allComb(k,1)).sessionNum;
        predCross_results.sessionNum_oblique_session     = data_info(allComb(k,2)).sessionNum;
        predCross_results.predResults                    = predResults;
        save(save_name, 'predCross_results')
    
    end
end
%% visualization for cross prediction


subjectCode = 'Ro';
dlist = dir(fullfile(save_folder_pred,sprintf('predCross_bestparams_monkey_%s_*',subjectCode)));
nSession = numel(dlist);

[fitCpredC_ll,fitOpredO_ll,...
 fitCpredO_ll,fitOpredC_ll,...
 fitCpredO_ll_flip,fitOpredC_ll_flip,...
 fitCpredO_ll_final,fitOpredC_ll_final] = deal(zeros(nSession,20));
for n = 1:nSession
    load(fullfile(save_folder_pred, dlist(n).name));

   
    fitCpredC_ll(n,:) = predCross_results.predResults.ll.fitCpredC';
    fitOpredO_ll(n,:) = predCross_results.predResults.ll.fitOpredO';


    fitCpredO_ll(n,:) = predCross_results.predResults.ll.fitCpredO';
    fitOpredC_ll(n,:) = predCross_results.predResults.ll.fitOpredC';


    fitCpredO_ll_flip(n,:) = predCross_results.predResults.ll.fitCpredO_flip';
    fitOpredC_ll_flip(n,:) = predCross_results.predResults.ll.fitOpredC_flip';

   
    % decide whether to flip the sign for this whole session
    if mean(fitCpredO_ll_flip(n,:)) > mean(fitCpredO_ll(n,:))
        fitCpredO_ll_final(n,:) = fitCpredO_ll_flip(n,:);
    else
        fitCpredO_ll_final(n,:)  = fitCpredO_ll(n,:);
    end

    if mean(fitOpredC_ll_flip(n,:)) > mean(fitOpredC_ll(n,:))
        fitOpredC_ll_final(n,:) = fitOpredC_ll_flip(n,:);
    else
        fitOpredC_ll_final(n,:)  = fitOpredC_ll(n,:);
    end
end
figure
subplot(1,2,1)
fig.plot_scatter_errorbar(fitCpredC_ll,fitOpredC_ll_final,...
            bpGlobal.color_list.color_cardinal,bpGlobal.color_list.color_cardinal_light);

[~,p,~,stats] = ttest(mean(fitCpredC_ll,2),mean(fitOpredC_ll_final,2));
title(sprintf(' $t(%d) = %.2f^{%s}$', stats.df,stats.tstat,fig.p2star(p)),...
 'Interpreter','latex','Color',bpGlobal.color_list.color_cardinal,'fontsize',14);
line([-1.1,-0.4],[-1.1,-0.4],'linewidth',3,'color','black','linestyle','--')
xlabel('LL (within-task)'); ylabel('LL (cross-task)')
set(gca,'FontSize',18)

subplot(1,2,2)
fig.plot_scatter_errorbar(fitOpredO_ll,fitCpredO_ll_final,...
            bpGlobal.color_list.color_oblique,bpGlobal.color_list.color_oblique_light);

[~,p,~,stats] = ttest(mean(fitOpredO_ll,2),mean(fitCpredO_ll_final,2));
title(sprintf(' $t(%d) = %.2f^{%s}$', stats.df,stats.tstat,fig.p2star(p)),...
 'Interpreter','latex','Color',bpGlobal.color_list.color_oblique,'fontsize',14);
line([-1.1,-0.4],[-1.1,-0.4],'linewidth',3,'color','black','linestyle','--')
xlabel('LL (within-task)'); ylabel('LL (cross-task)')
set(gca,'FontSize',18)

sgtitle(sprintf('Using best parameters for %s',subjectCode),'fontsize',20,'fontweight','bold')