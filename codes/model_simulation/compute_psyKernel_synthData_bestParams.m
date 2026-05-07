clear all
clc
close all
%%%% By Shizhao Liu 03/07/2026
%%%%% Cpmputr psychometric kernels of synthetic data generated with the
%%%%% best few sets of parameters.
%%%%% Also runs cross-prediction
global bpGlobal
bpGratingFCGlobal()
%% Compute psychometric kernel
dataFolder_home = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psyKernel';
doKernel = 0;
if doKernel
    save_folder = '../../results/model_behavior/psyKernel_bestParams';
    nSet = 10; % 10 sets of parameters for each animal
    
    animal_list = {'monkeyR';'monkeyG'};
    energy_use = 'proj';
    for i = 1:numel(animal_list)
        for j = 1:nSet
            
            dataFolder = fullfile(dataFolder_home,sprintf('best_parameters_%s_set_%d_contrast6', animal_list{i},j));
    
            data_list  = dir(fullfile(dataFolder,'*.mat'));
            
            psyKernel_model = struct();
            save_name       = fullfile(save_folder,sprintf('psyKernels_best_parameters_%s_set_%d', animal_list{i},j));
            for n = 1:numel(data_list)
                fprintf('Running animal %d, para set %d, session %d \n', i,j,n)
                % load
                load(fullfile(dataFolder, data_list(n).name));
    
                tokens = regexp(data_list(n).name, 'Single_image_task_([a-zA-Z0-9]+)_delta_([-]?[\d_]+)_prior_([-]?[\d_]+)_session_([\d_]+)', 'tokens');
                % Convert extracted strings back to numbers
                extracted_params        = tokens{1}; % Extract matched tokens
               
                task_str                = extracted_params{1};
                delta_str               = strrep(extracted_params{2}, '_', '.');
                prior_str               = strrep(extracted_params{3}, '_', '.');
                session_num_str         = extracted_params{4};
    
                % compute
                [map_params,orientationsDEG] = util_it.compute_psyKernel_model(data_use, energy_use);
                % save
                psyKernel_model(n).subjectCode      = animal_list{i};
                psyKernel_model(n).task             = task_str;
                psyKernel_model(n).delta            = str2double(delta_str);
                psyKernel_model(n).prior            = str2double(prior_str);    
                psyKernel_model(n).repeat           = str2double(session_num_str);
                psyKernel_model(n).w_ori            = map_params.w_ori;
                psyKernel_model(n).w_time           = map_params.w_time;
                psyKernel_model(n).orientationsDEG  = orientationsDEG;
            
                psyKernel_model(n).energy_use       = energy_use;
            end
    
        save(save_name, 'psyKernel_model');
           
            
    
        end
    end
end

%% check results
%%%%%% Compute correlations between kernels
%%%%% load empirical kernel

psyKernel_Ro = load('/Users/liushizhao/projects_local/bpGratingEx/results/behavior/Rolo_psyKernel_table_final');
psyKernel_Gr = load('/Users/liushizhao/projects_local/bpGratingEx/results/behavior/Gremlin_psyKernel_table');

sessionName_list  = bpGlobal.rolo.session_list.switching;
idx_list          = ismember({psyKernel_Ro.psyKernel_table(:).sessionName},sessionName_list);

tmp                         = {psyKernel_Ro.psyKernel_table(idx_list).w_ori_cardinal};
w_ori_cardinal_Ro = cat(2,tmp{:});
tmp                         = {psyKernel_Ro.psyKernel_table(idx_list).w_ori_oblique};
w_ori_oblique_Ro  = cat(2,tmp{:});

sessionName_list  = bpGlobal.gremlin.session_list.interleaved_real;
idx_list          = ismember({psyKernel_Gr.psyKernel_table(:).sessionName},sessionName_list);

tmp                         = {psyKernel_Gr.psyKernel_table(idx_list).w_ori_cardinal};
w_ori_cardinal_Gr = cat(2,tmp{:});
tmp                         = {psyKernel_Gr.psyKernel_table(idx_list).w_ori_oblique};
w_ori_oblique_Gr  = cat(2,tmp{:});


orientationsDEG                  = [0:15:165];
w_ori_cardinal_ideal            = 0.4*sin(pi * (orientationsDEG - 45)/90)';
w_ori_oblique_ideal             = 0.4*sin(pi * (orientationsDEG - 90)/90)';


subject_list = {'monkeyR';'monkeyG'};
data_folder  = '../../results/model_behavior/psyKernel_bestParams';

k  = 1;
for i = 1:2
    for i_set = 1:10
        subjectCode = subject_list{i};
   
        data_name = fullfile(data_folder, sprintf('psyKernels_best_parameters_%s_set_%d', subjectCode, i_set));
        load(data_name);
        
        
        idx_cardinal = strcmp({psyKernel_model(:).task},'cardinal');
        idx_oblique = strcmp({psyKernel_model(:).task},'oblique');
        prior_cardinal  =  psyKernel_model(find(idx_cardinal,1)).prior;
        prior_oblique   =  psyKernel_model(find(idx_oblique,1)).prior;
        delta_cardinal  =  psyKernel_model(find(idx_cardinal,1)).delta;
        delta_oblique   = psyKernel_model(find(idx_oblique,1)).delta;
        tmp = {psyKernel_model(idx_cardinal).w_ori};
        w_ori_cardinal_model = cat(2, tmp{:});
        tmp = {psyKernel_model(idx_oblique).w_ori};
        w_ori_oblique_model = cat(2, tmp{:});
        
        switch subjectCode
            case 'monkeyR'
                w_ori_cardinal_empirical = w_ori_cardinal_Ro;
                w_ori_oblique_empirical = w_ori_oblique_Ro;
            case 'monkeyG'
                w_ori_cardinal_empirical = w_ori_cardinal_Gr;
                w_ori_oblique_empirical = w_ori_oblique_Gr;
        end
        %%% corr_model_empirical
        [r_model_empirical_cardinal_avg, r_model_empirical_cardinal_combinations] = ...
                compute_corr_kernels(w_ori_cardinal_model, w_ori_cardinal_empirical);
        [r_model_empirical_oblique_avg, r_model_empirical_oblique_combinations] = ...
                compute_corr_kernels(w_ori_oblique_model, w_ori_oblique_empirical);
        
        %%% corr_model_ideal
        [r_model_ideal_cardinal_avg, r_model_ideal_cardinal_combinations] = ...
                compute_corr_kernels(w_ori_cardinal_model, w_ori_cardinal_ideal);
        [r_model_ideal_oblique_avg, r_model_ideal_oblique_combinations] = ...
                compute_corr_kernels(w_ori_oblique_model, w_ori_oblique_ideal);

        %%% corr_empirical_ideal
        [r_empirical_ideal_cardinal_avg, ~] = ...
                compute_corr_kernels(w_ori_cardinal_empirical, w_ori_cardinal_ideal);
        [r_empirical_ideal_oblique_avg, ~] = ...
                compute_corr_kernels(w_ori_oblique_empirical, w_ori_oblique_ideal);

        %%%% corr_model_ideal_cross
        [r_model_ideal_cross_cardinal_avg, ~] = ...
                compute_corr_kernels(w_ori_cardinal_model, w_ori_oblique_ideal);
        [r_model_ideal_cross_oblique_avg, ~] = ...
                compute_corr_kernels(w_ori_oblique_model, w_ori_cardinal_ideal);

        %%%% corr_empirical_ideal_cross
        [r_empirical_ideal_cross_cardinal_avg, ~] = ...
                compute_corr_kernels(w_ori_cardinal_empirical, w_ori_oblique_ideal);
        [r_empirical_ideal_cross_oblique_avg, ~] = ...
                compute_corr_kernels(w_ori_oblique_empirical, w_ori_cardinal_ideal);

        
        psyKernel_corr(k).subjectCode                               = subjectCode;
        psyKernel_corr(k).i_set                                     = i_set;
        psyKernel_corr(k).prior_cardinal                            = prior_cardinal;
        psyKernel_corr(k).prior_oblique                             = prior_oblique;
        psyKernel_corr(k).delta_cardinal                            = delta_cardinal;
        psyKernel_corr(k).delta_oblique                             = delta_oblique;
        psyKernel_corr(k).r_model_empirical_cardinal_avg            = r_model_empirical_cardinal_avg;
        psyKernel_corr(k).r_model_empirical_cardinal_combinations   = r_model_empirical_cardinal_combinations;
        
        psyKernel_corr(k).r_model_empirical_oblique_avg             = r_model_empirical_oblique_avg;
        psyKernel_corr(k).r_model_empirical_oblique_combinations    = r_model_empirical_oblique_combinations;
        
        psyKernel_corr(k).r_model_ideal_cardinal_avg                = r_model_ideal_cardinal_avg;
        psyKernel_corr(k).r_model_ideal_cardinal_combinations       = r_model_ideal_cardinal_combinations;
        psyKernel_corr(k).r_model_ideal_oblique_avg                 = r_model_ideal_oblique_avg;
        psyKernel_corr(k).r_model_ideal_oblique_combinations        = r_model_ideal_oblique_combinations;

        psyKernel_corr(k).r_empirical_ideal_cardinal_avg            = r_empirical_ideal_cardinal_avg;
        psyKernel_corr(k).r_empirical_ideal_oblique_avg             = r_empirical_ideal_oblique_avg;


        psyKernel_corr(k).r_model_ideal_cross_cardinal_avg          = r_model_ideal_cross_cardinal_avg;
        psyKernel_corr(k).r_model_ideal_cross_oblique_avg           = r_model_ideal_cross_oblique_avg;
        psyKernel_corr(k).r_empirical_ideal_cross_cardinal_avg      = r_empirical_ideal_cross_cardinal_avg;
        psyKernel_corr(k).r_empirical_ideal_cross_oblique_avg       = r_empirical_ideal_cross_oblique_avg;


        k = k+1;
    end
end
%% visualization of correlation values
%%%%% scatter, corr_model_empirical_cardinal, corr_model_empirical_oblique
%%%%% scatter, corr_model_ideal_cardinal, corr_model_ideal_oblique. With
%%%%% the corr_empirical_ideal marked
subjectCode = 'monkeyG';

figure
subplot(1,2,1); hold on

idx = strcmp({psyKernel_corr(:).subjectCode}, subjectCode);
x_avg = [psyKernel_corr(idx).r_model_empirical_cardinal_avg];
y_avg = [psyKernel_corr(idx).r_model_empirical_oblique_avg];
i_set = [psyKernel_corr(idx).i_set];
for n = 1:numel(x_avg)
    text(x_avg(n), y_avg(n), sprintf('%d', i_set(n)),'fontsize',16)
end
set(gca,'fontsize',16)
xlabel('Cardinal');
ylabel('Oblique')
title([subjectCode, ' corr-model-empirical'])
addCircularGrid([0,0], [0:0.05:1.2]);
xlim([0.5,1]); ylim([0.5,1])

subplot(1,2,2); hold on

idx = strcmp({psyKernel_corr(:).subjectCode}, subjectCode);
x_avg = [psyKernel_corr(idx).r_model_ideal_cardinal_avg];
y_avg = [psyKernel_corr(idx).r_model_ideal_oblique_avg];
i_set = [psyKernel_corr(idx).i_set];
scatter(x_avg, y_avg)
for n = 1:numel(x_avg)
    text(x_avg(n), y_avg(n), sprintf('%d', i_set(n)),'fontsize',16)
end
set(gca,'fontsize',16)
xlabel('Cardinal');
ylabel('Oblique')
title([subjectCode, ' corr-model-ideal'])

r_empirical_ideal_cardinal_avg = psyKernel_corr(find(idx,1)).r_empirical_ideal_cardinal_avg;
r_empirical_ideal_oblique_avg = psyKernel_corr(find(idx,1)).r_empirical_ideal_oblique_avg;

plot(r_empirical_ideal_cardinal_avg, r_empirical_ideal_oblique_avg, '*', 'MarkerSize',20,'Color','red');

addCircularGrid([r_empirical_ideal_cardinal_avg,r_empirical_ideal_oblique_avg], [0:0.05:1]);
xlim([0.5,1]); ylim([0.5,1])

%% visualization of the model kernels
kernel_type = 'spatial';
subjectCode = 'monkeyG';
i_set = 8;
plot_psykernel(subjectCode, i_set, kernel_type)
%% Runs cross prediction
doPred = 1;
if doPred
    monkeyR_set_list = [7,8];
    monkeyG_set_list = [4,8];
  
    nFold = 20;
    nCombination_Max = 25;
    dataFolder_home = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psyKernel';
    animal_list = {'monkeyR';'monkeyG'};
    
    for i = 1:numel(animal_list)
        subjectCode = animal_list{i};
        switch subjectCode
            case 'monkeyR'
                set_list = monkeyR_set_list;
            case 'monkeyG'
                set_list = monkeyG_set_list;
        end
        nSet = numel(set_list);
        for j = 1:nSet
            i_set = set_list(j);
            data_folder = fullfile(dataFolder_home,sprintf('best_parameters_%s_set_%d_contrast6', subjectCode,i_set));
            save_folder_pred = fullfile('../../results/model_behavior/predCross', sprintf('bestparams_%s_set_%d',subjectCode,i_set));
            if ~isfolder(save_folder_pred)
                mkdir(save_folder_pred)
            end
            data_list  = dir(fullfile(data_folder,'*.mat'));

            for n = 1:numel(data_list)
          
                tokens = regexp(data_list(n).name, '([a-zA-Z0-9]+)_image_task_([a-zA-Z0-9]+)_delta_([-]?[\d_]+)_prior_([\d_]+)_session_([\d_]+)', 'tokens');
                % Convert extracted strings back to numbers
                extracted_params        = tokens{1}; % Extract matched tokens
                task_mode_str           = extracted_params{1};
                imagetask_str           = extracted_params{2};
                delta_str               = strrep(extracted_params{3}, '_', '.');
                prior_str               = strrep(extracted_params{4}, '_', '.');
                session_num_str         = extracted_params{5};
            
            
            
                data_info(n).filename           = data_list(n).name;
                data_info(n).mode               = task_mode_str;
                data_info(n).image_task         = imagetask_str;
                data_info(n).delta              = str2double(delta_str);
                data_info(n).prior              = str2double(prior_str);
                data_info(n).sessionNum         = str2double(session_num_str);
            end

            idx_cardinal = find(strcmp({data_info(:).image_task},'cardinal'));
            idx_oblique  = find(strcmp({data_info(:).image_task},'oblique'));
          
    
            [A, B] = ndgrid(idx_cardinal, idx_oblique);
            allComb = [A(:) B(:)];
    
            nCombination = min([size(allComb,1), nCombination_Max]);
            for k = 1:nCombination
                save_name           = fullfile(save_folder_pred,sprintf('predCross_bestparams_monkey_%s_combination_%d.mat',...
                    subjectCode, k));
                if isfile(save_name)
                    continue
                end
                fprintf('Running monkey %d, set %d, combination %d/%d \n',i,j,k, nCombination) ;
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
    end
end
%% visualization of the cross-prediction results
subjectCode = 'monkeyG';
i_set = 8;
plot_predCross_results(subjectCode, i_set);

%% check how good these selected sets are in terms of matching neural data

load('../../results/model_match_neural_rolo_moreRepeats.mat');
load('../../results/model_match_neural_gremlin_moreRepeats.mat');

[~, index_rolo] = sort(neural_match_result_rolo.dist_feature, 'ascend');
[~, index_gremlin] = sort(neural_match_result_gremlin.dist_feature, 'ascend');

monkeyR_set_list = [7,8];
monkeyG_set_list = [4,8];
nPlot = 30;

figure;
subplot(2,3,[1,2]); hold on
plot(neural_match_result_rolo.dist_feature(index_rolo),'LineWidth',2);
set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); xlim([0,900])
title('Monkey R')
subplot(2,3,3); hold on
plot(neural_match_result_rolo.dist_feature(index_rolo(1:nPlot)),'-o');
for n = 1:numel(monkeyR_set_list)
    i = monkeyR_set_list(n);
    plot(i, neural_match_result_rolo.dist_feature(index_rolo(i)),'*','color','black');
end
set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); 


subplot(2,3,[4,5]); hold on
plot(neural_match_result_gremlin.dist_feature(index_gremlin),'LineWidth',2);
set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); xlim([0,900])
title('Monkey G')
subplot(2,3,6); hold on
plot(neural_match_result_gremlin.dist_feature(index_gremlin(1:nPlot)),'-o');
for n = 1:numel(monkeyR_set_list)
    i = monkeyR_set_list(n);
    plot(i, neural_match_result_gremlin.dist_feature(index_gremlin(i)),'*','color','black');
end
set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); 
%% helper functions
function [r_avgs, r_combinations] = compute_corr_kernels(w_ori_1, w_ori_2)
if size(w_ori_1,1) ~= size(w_ori_2,1)
    error('Size of psychometric kernels need to match')
end
w_ori_1_avg = mean(w_ori_1, 2);
w_ori_2_avg = mean(w_ori_2, 2);

r_avgs = corr(w_ori_1_avg, w_ori_2_avg);

N_1 = size(w_ori_1, 2);
N_2 = size(w_ori_2, 2);
r_combinations = zeros(N_1, N_2);
for n1 = 1:N_1
    for n2 = 1:N_2
        r_combinations(n1, n2) = corr(w_ori_1(:, n1), w_ori_2(:, n2)); 
    end
end
r_combinations = r_combinations(:);

end

function plot_predCross_results(subjectCode, i_set)
    global bpGlobal
    home_folder = '../../results/model_behavior/predCross';
    data_folder = fullfile(home_folder, sprintf('bestparams_%s_set_%d',subjectCode, i_set));
    dlist   = dir(fullfile(data_folder,'*.mat'));
    nSession = numel(dlist);
    
    [fitCpredC_ll,fitOpredO_ll,...
     fitCpredO_ll,fitOpredC_ll,...
     fitCpredO_ll_flip,fitOpredC_ll_flip,...
     fitCpredO_ll_final,fitOpredC_ll_final] = deal(zeros(nSession,20));
    for n = 1:nSession
        load(fullfile(dlist(n).folder, dlist(n).name));
    
       
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

    sgtitle(sprintf('Parameters set %d for %s',i_set, subjectCode),'fontsize',20,'fontweight','bold')

end

function plot_psykernel(subjectCode, i_set, kernel_type)
data_folder = '../../results/model_behavior/psyKernel_bestParams';
data_name = fullfile(data_folder, sprintf('psyKernels_best_parameters_%s_set_%d', subjectCode, i_set));
load(data_name);
for n = 1:numel(psyKernel_model)
    %%%% add these fields just for compaibility of the plot function
    psyKernel_model(n).mode = 'Single';
    psyKernel_model(n).image_task = psyKernel_model(n).task;
end
idx_cardinal    = strcmp({psyKernel_model(:).task},'cardinal');
idx_oblique     = strcmp({psyKernel_model(:).task},'oblique');

prior_cardinal_plot = psyKernel_model(find(idx_cardinal,1)).prior;
delta_cardinal_plot = psyKernel_model(find(idx_cardinal,1)).delta;

prior_oblique_plot = psyKernel_model(find(idx_oblique,1)).prior;
delta_oblique_plot = psyKernel_model(find(idx_oblique,1)).delta;


figure
subplot(1,2,1)
plotOptions.mode            = 'Single';
plotOptions.image_task      = 'cardinal';
plotOptions.prior_task      = prior_cardinal_plot;
plotOptions.delta           = delta_cardinal_plot;
plotOptions.plotIndividual  = true;
fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, plotOptions);
title(sprintf('Delta = %.2f,  Prior = %.1f', ...
    plotOptions.delta, plotOptions.prior_task));

subplot(1,2,2)
plotOptions.mode            = 'Single';
plotOptions.image_task      = 'oblique';
plotOptions.prior_task      = prior_oblique_plot;
plotOptions.delta           = delta_oblique_plot;
plotOptions.plotIndividual  = true;
fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, plotOptions);
title(sprintf('Delta = %.2f,  Prior = %.1f', ...
    plotOptions.delta, plotOptions.prior_task));

sgtitle(sprintf('Parameters set %d for %s',i_set, subjectCode),'fontsize',20,'fontweight','bold')
end



function addCircularGrid(center, radii)

x0 = center(1);
y0 = center(2);

theta = linspace(0,2*pi,300);

hold on

for r = radii
    x = x0 + r*cos(theta);
    y = y0 + r*sin(theta);
    plot(x,y,'Color','k','LineStyle','--','LineWidth',1)
end

end