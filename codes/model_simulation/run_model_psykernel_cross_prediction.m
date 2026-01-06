clear all
clc
close all
%%
home_folder     = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel';

version_name    = 'moreRepeats_contrast_6';
data_folder     = fullfile(home_folder, version_name);

save_folder     = fullfile('../../results/model_behavior/predCross', version_name);
if ~isfolder(save_folder)
    mkdir(save_folder)
end


%% compute: prediction and cross prediction
doCompute = 1;
if doCompute
    nCombination = 30;
    nFold = 20;
    data_list = dir(fullfile(data_folder,'*.mat'));
    
        
    data_info = struct();
    for n = 1:numel(data_list)
      
        tokens = regexp(data_list(n).name, '([a-zA-Z0-9]+)_image_task_([a-zA-Z0-9]+)_prior_cardinal_([-]?[\d_]+)_prior_oblique_([\d_]+)_session_([\d_]+)', 'tokens');
        % Convert extracted strings back to numbers
        extracted_params        = tokens{1}; % Extract matched tokens
        task_mode_str           = extracted_params{1};
        imagetask_str           = extracted_params{2};
        prior_cardinal_str      = strrep(extracted_params{3}, '_', '.');
        prior_oblique_str       = strrep(extracted_params{4}, '_', '.');
        session_num_str         = extracted_params{5};
    
    
        data_info(n).filename         =  data_list(n).name;
        data_info(n).mode             = task_mode_str;
        data_info(n).image_task       = imagetask_str;
        data_info(n).prior_cardinal   = str2double(prior_cardinal_str);
        data_info(n).prior_oblique    = str2double(prior_oblique_str);
        data_info(n).sessionNum       = str2double(session_num_str);
    end
    %%% find all unique combination of priors
    
    idx_cardinal = strcmp({data_info(:).image_task},'cardinal');
    prior_cardinal_cardinal = [data_info(idx_cardinal).prior_cardinal]';
    prior_cardinal_oblique  = [data_info(idx_cardinal).prior_oblique]';
    unique_prior_pairs_cardinal = unique([prior_cardinal_cardinal prior_cardinal_oblique], 'rows');
    
    %unique_prior_pairs_cardinal = [1,0]; % over write this
    
    
    idx_oblique = strcmp({data_info(:).image_task},'oblique');
    prior_oblique_cardinal = [data_info(idx_oblique).prior_cardinal]';
    prior_oblique_oblique  = [data_info(idx_oblique).prior_oblique]';
    unique_prior_pairs_oblique = unique([prior_oblique_cardinal prior_oblique_oblique], 'rows');
    
    
    nCardinal_prior = size(unique_prior_pairs_cardinal,1);
    nOblique_prior  = size(unique_prior_pairs_oblique, 1); 
    for i = 1:nCardinal_prior
        for j = 1:nOblique_prior
            %%%% the first cardinal refers to cardinal task, the second refers to
            %%%% the prior set-up
            prior_cardinal_cardinal = unique_prior_pairs_cardinal(i,1);
            prior_cardinal_oblique  = unique_prior_pairs_cardinal(i,2);
            prior_oblique_cardinal  = unique_prior_pairs_oblique(j,1);
            prior_oblique_oblique   = unique_prior_pairs_oblique(j,2);
    
            idx_cardinal = find(strcmp({data_info(:).image_task},'cardinal') & ...
                        [data_info(:).prior_cardinal] == prior_cardinal_cardinal & ...
                        [data_info(:).prior_oblique] == prior_cardinal_oblique);
            idx_oblique = find(strcmp({data_info(:).image_task},'oblique') & ...
                        [data_info(:).prior_cardinal] == prior_oblique_cardinal & ...
                        [data_info(:).prior_oblique] == prior_oblique_oblique);
            
            [A, B] = ndgrid(idx_cardinal, idx_oblique);
            allComb = [A(:) B(:)];
            nPick = min(nCombination, size(allComb,1));   % safety
            sel = randperm(size(allComb,1), nPick);
            randComb = allComb(sel,:);
    
            
           
        
            prior_cardinal_cardinal_str  = strrep(sprintf('%.1f',prior_cardinal_cardinal), '.', '_');
            prior_cardinal_oblique_str   = strrep(sprintf('%.1f',prior_cardinal_oblique), '.', '_');
            prior_oblique_cardinal_str   = strrep(sprintf('%.1f',prior_oblique_cardinal), '.', '_');
            prior_oblique_oblique_str    = strrep(sprintf('%.1f',prior_oblique_oblique), '.', '_');
    
            
    
            for k = 20:size(randComb, 1)
                save_name           = fullfile(save_folder,sprintf('predCross_prior_cardinal_session_%s_%s_prior_oblique_session_%s_%s_combination_%d.mat',...
                                            prior_cardinal_cardinal_str, prior_cardinal_oblique_str, ...
                                            prior_oblique_cardinal_str,  prior_oblique_oblique_str, k)); 
                if isfile(save_name)
                    continue
                end
                fprintf('Running cardinal %d/%d, oblique %d/%d, combination %d/%d \n',...
                        i, nCardinal_prior, j, nOblique_prior, k, nCombination );
                file_name_cardinal  = fullfile(data_folder, data_info(randComb(k,1)).filename);
                file_name_oblique   = fullfile(data_folder, data_info(randComb(k,2)).filename);
                
                dat_cardinal        = load(file_name_cardinal);
                dat_oblique         = load(file_name_oblique);
                predResults         = run_cross_predict_model(dat_cardinal,dat_oblique, nFold);
                
                predCross_results   = struct();
                predCross_results.mode_cardinal_session          = data_info(randComb(k,1)).mode;
                predCross_results.mode_oblique_session           = data_info(randComb(k,2)).mode;
                predCross_results.prior_cardinal_session         = unique_prior_pairs_cardinal(i,:);
                predCross_results.prior_oblique_session          = unique_prior_pairs_oblique(j,:);
                predCross_results.sessionNum_cardinal_session    = data_info(randComb(k,1)).sessionNum;
                predCross_results.sessionNum_oblique_session     = data_info(randComb(k,2)).sessionNum;
                predCross_results.predResults                    = predResults;
                save(save_name, 'predCross_results')
    
            end
            
        end
    end
end
%% visualization
dlist = dir(fullfile(save_folder,'*.mat'));
nSession = numel(dlist);

[fitCpredC_ll,fitOpredO_ll,...
 fitCpredO_ll,fitOpredC_ll,...
 fitCpredO_ll_flip,fitOpredC_ll_flip,...
 fitCpredO_ll_final,fitOpredC_ll_final] = deal(zeros(nSession,20));

[prior_cardinal_session, prior_oblique_session] = deal(zeros(nSession,2));
for n = 1:nSession
    load(fullfile(save_folder, dlist(n).name));

    prior_cardinal_session(n,:) = predCross_results.prior_cardinal_session;
    prior_oblique_session(n,:) = predCross_results.prior_oblique_session;
    
    fitCpredC_ll(n,:) = predCross_results.predResults.ll.fitCpredC';
    fitOpredO_ll(n,:) = predCross_results.predResults.ll.fitOpredO';


    fitCpredO_ll(n,:) = predCross_results.predResults.ll.fitCpredO';
    fitOpredC_ll(n,:) = predCross_results.predResults.ll.fitOpredC';


    fitCpredO_ll_flip(n,:) = predCross_results.predResults.ll.fitCpredO_flip';
    fitOpredC_ll_flip(n,:) = predCross_results.predResults.ll.fitOpredC_flip';

    fitCpredO_ll_final(n,:) = fitCpredO_ll(n,:);
    fitOpredC_ll_final(n,:) = fitOpredC_ll(n,:);
    % % decide whether to flip the sign for this whole session
    % if mean(fitCpredO_ll_flip(n,:)) > mean(fitCpredO_ll(n,:))
    %     fitCpredO_ll_final(n,:) = fitCpredO_ll_flip(n,:);
    % else
    %     fitCpredO_ll_final(n,:)  = fitCpredO_ll(n,:);
    % end
    % 
    % if mean(fitOpredC_ll_flip(n,:)) > mean(fitOpredC_ll(n,:))
    %     fitOpredC_ll_final(n,:) = fitOpredC_ll_flip(n,:);
    % else
    %     fitOpredC_ll_final(n,:)  = fitOpredC_ll(n,:);
    % end
end
%%
global bpGlobal
bpGratingFCGlobal();

unique_prior_pairs_cardinal = unique(prior_cardinal_session, 'rows');
unique_prior_pairs_oblique  = unique(prior_oblique_session, 'rows');

% [~,idx_sort] = sort(unique_prior_pairs_cardinal(:,1),'descend');
% unique_prior_pairs_cardinal = unique_prior_pairs_cardinal(idx_sort,:);
nCardinal_prior = size(unique_prior_pairs_cardinal,1);
nOblique_prior  = size(unique_prior_pairs_oblique, 1); 
for i = 1:nCardinal_prior
    for j = 1:nOblique_prior
        idx = ismember(prior_cardinal_session, unique_prior_pairs_cardinal(i,:), 'rows') & ...
               ismember(prior_oblique_session, unique_prior_pairs_oblique(j,:), 'rows');

        fitCpredC_ll_plot = fitCpredC_ll(idx,:);
        fitOpredC_ll_plot = fitOpredC_ll_final(idx,:);
    
        fitOpredO_ll_plot = fitOpredO_ll(idx,:);
        fitCpredO_ll_plot = fitCpredO_ll_final(idx,:);
        
        subplot(4,8,(i-1) * 8 + (j-1) * 2 + 1)
        fig.plot_scatter_errorbar(fitCpredC_ll_plot,fitOpredC_ll_plot,...
            bpGlobal.color_list.color_cardinal,bpGlobal.color_list.color_cardinal_light);
        [~,p,~,stats] = ttest(mean(fitCpredC_ll_plot,2),mean(fitOpredC_ll_plot,2));
        title(sprintf('prior = %.1f \n $t(%d) = %.2f^{%s}$',unique_prior_pairs_cardinal(i,1), stats.df,stats.tstat,fig.p2star(p)),...
         'Interpreter','latex','Color',bpGlobal.color_list.color_cardinal,'fontsize',14);


        line([-1.1,-0.4],[-1.1,-0.4],'linewidth',3,'color','black','linestyle','--')
        xlabel('LL (within-task)'); ylabel('LL (cross-task)')

        subplot(4,8,(i-1) * 8 + j * 2)
        fig.plot_scatter_errorbar(fitOpredO_ll_plot,fitCpredO_ll_plot,...
            bpGlobal.color_list.color_oblique,bpGlobal.color_list.color_oblique_light);
        [~,p,~,stats] = ttest(mean(fitOpredO_ll_plot,2),mean(fitCpredO_ll_plot,2));
        title(sprintf('prior = %.1f \n $t(%d) = %.2f^{%s}$',unique_prior_pairs_oblique(j,2),stats.df,stats.tstat,fig.p2star(p)),...
         'Interpreter','latex','Color',bpGlobal.color_list.color_oblique,'fontsize',14);

        line([-1.1,-0.4],[-1.1,-0.4],'linewidth',3,'color','black','linestyle','--')
        xlabel('LL (within-task)'); ylabel('LL (cross-task)')

    end
end
%%

function predResults = run_cross_predict_model(dat_cardinal,dat_oblique, nFold)

energy_use = 'proj';
kernelOptions.orientationsDEG = [0:15:165]; % orientation bin: 0 to 180, 15 degree increment
kernelOptions.filterMode ='hard';
kernelOptions.filterSDDeg = 15;
fittingKernelParams.hypers  = [0.5 6 60 0.15 0.06];
fittingKernelParams.nLapse  = 2;
[X_cardinal_norm, choice_cardinal, contrast_cardinal] = util_it.extract_X_psyKernel(dat_cardinal.data_use, energy_use);
[X_oblique_norm,  choice_oblique, contrast_oblique] = util_it.extract_X_psyKernel(dat_oblique.data_use, energy_use);

X.cardinal            = X_cardinal_norm;
X.oblique             = X_oblique_norm;
choice.cardinal       = choice_cardinal;
choice.oblique        = choice_oblique;
choice_list.cardinal  = unique(choice_cardinal);
choice_list.oblique   = unique(choice_oblique);
signal_level.cardinal = contrast_cardinal;
signal_level.oblique  = contrast_oblique;

predResults = psyKernel.fitSPTP_kernel_crosstest(X,choice,kernelOptions,choice_list,signal_level,nFold,fittingKernelParams);

end