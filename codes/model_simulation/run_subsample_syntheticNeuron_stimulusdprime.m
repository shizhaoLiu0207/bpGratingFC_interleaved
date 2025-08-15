clear all
clc
close all
%%%% This script tests subsampling model neurons so that the distribution
%%%% of their stimulus dprime matches empirical data

%% get neural response of synthetic neurons and their stimulus dprime
%%%%% These lists of stimulus_contrast, b_PF and image_task are consistent
%%%%% with what's used to generate large-scale synthetic data
b_PF_list           = [0.8, 0];
nNeuron             = 128;

stimulus_signal    = 12;

nTrial = 1000;
k = 1;
projective_fields = cell(numel(b_PF_list),1);

for i = 1:numel(b_PF_list)
    %for n = 1:numel(stimulus_signal_list)
        % i = 1;
        % n = 1;
        %%%%%% Fake task parameters. Fixed. It does not affect neuron's feedforward
        %%%%%% stimulus tuning. This is just to generate P and projective_fields
        image_task              = 'cardinal';
        prior_task              = [0.5,0.5];
        delta                   = 0.08;
        stimulus_contrast       = [0,0];
        b_sample                    = b_PF_list(i);
        P = S_Exp_Para('run-interleaved-nonuniform', 'I.stimulus_contrast',stimulus_contrast,...
            'G.prior_task',prior_task,...
            'I.image_task',image_task,...
            'G.b_PF',b_sample,...
            'G.delta',delta,...
            'G.dimension_X',nNeuron);
        projective_fields{i} = C_Projection(P.G.fct, P.G.nx, P.G.dimension_X, P.G.dimension_G, P.G.number_locations, P.G.b_PF);

        im_type     = P.I.fct;
        n_locs      = P.G.number_locations;
        im_height   = P.G.ny;
       % nNeuron     = size(projective_fields{i}.G,2);
        
        %%%%%%%% Real experiment. Generate neuron response and dprime
        %%%% generate parameter struct



        %stimulus_signal     = stimulus_signal_list(n);


        [response_cardinal_1,response_cardinal_2,response_oblique_1,response_oblique_2] = deal(zeros(nTrial,nNeuron));
        for t = 1:nTrial
            %%%%% cardinal image
            image_task  = 'cardinal';
            stimulus_contrast = [stimulus_signal,0];
            stim_cardinal_1 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);
            image_task  = 'cardinal';
            stimulus_contrast = [0,stimulus_signal];
            stim_cardinal_2 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);

            %%%%% cardinal signal
            response_cardinal_1(t,:) = arrayfun(@(n)stim_cardinal_1(:)' * projective_fields{i}.G(:,n), [1:nNeuron]);
            response_cardinal_2(t,:) = arrayfun(@(n)stim_cardinal_2(:)' * projective_fields{i}.G(:,n), [1:nNeuron]);
            %cardinal_signal(i,t) = mean(abs(response_cardinal_1 - response_cardinal_2));


            %%%%% oblique image
            image_task  = 'oblique';
            stimulus_contrast = [stimulus_signal,0];
            stim_oblique_1 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);
            stimulus_contrast = [0,stimulus_signal];
            stim_oblique_2 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);
            %%%%% oblique signal
            response_oblique_1(t,:) = arrayfun(@(n)stim_oblique_1(:)' * projective_fields{i}.G(:,n), [1:nNeuron]);
            response_oblique_2(t,:) = arrayfun(@(n)stim_oblique_2(:)' * projective_fields{i}.G(:,n), [1:nNeuron]);
            %oblique_signal(i,t) = mean(abs(response_oblique_1 - response_oblique_2));




        end
        dprime_cardinal =  abs(mean(response_cardinal_1,1) - mean(response_cardinal_2,1)) ./ ...
            sqrt((var(response_cardinal_1, [], 1) + var(response_cardinal_2, [], 1)) / 2);
        dprime_oblique = abs(mean(response_oblique_1,1) - mean(response_oblique_2,1)) ./ ...
            sqrt((var(response_oblique_1, [], 1) + var(response_oblique_2, [], 1)) / 2);

        dprime_cardinal_sign =  (mean(response_cardinal_1,1) - mean(response_cardinal_2,1)) ./ ...
            sqrt((var(response_cardinal_1, [], 1) + var(response_cardinal_2, [], 1)) / 2);
        dprime_oblique_sign = (mean(response_oblique_1,1) - mean(response_oblique_2,1)) ./ ...
            sqrt((var(response_oblique_1, [], 1) + var(response_oblique_2, [], 1)) / 2);

        results_dprime(k).b_PF                      = b_sample;
        results_dprime(k).projective_fields         = projective_fields{i};
        results_dprime(k).stimulus_signal           = stimulus_signal;
        results_dprime(k).dprime_cardinal           = dprime_cardinal;
        results_dprime(k).dprime_oblique            = dprime_oblique;
        results_dprime(k).dprime_cardinal_sign      = dprime_cardinal_sign;
        results_dprime(k).dprime_oblique_sign       = dprime_oblique_sign;
        % results_dprime(k).dprime_cardinal_norm      = dprime_cardinal / (2 * stimulus_signal);
        % results_dprime(k).dprime_oblique_norm       = dprime_oblique / (2 * stimulus_signal);
        % results_dprime(k).dprime_cardinal_sign_norm = dprime_cardinal_sign / (2 * stimulus_signal);
        % results_dprime(k).dprime_oblique_sign_norm  = dprime_oblique_sign / (2 * stimulus_signal);
        k = k +1;

    %end
end

%% dprimes from empirical data and tuning
load('../../results/neural/stimulus_dprime_twoanimals_interleaved_all_trials_coef1_hVis2_FR1.mat');

idx_rolo_cardinal = contains({results_fprime_choice_all(:).sessionStr},'Ro') & ...
                strcmp({results_fprime_choice_all(:).task},'cardinal') & ...
                cell2mat({results_fprime_choice_all(:).cohr_level}) == 10;

idx_rolo_oblique = contains({results_fprime_choice_all(:).sessionStr},'Ro') & ...
                strcmp({results_fprime_choice_all(:).task},'oblique') & ...
                cell2mat({results_fprime_choice_all(:).cohr_level}) == 10;

idx_gremlin_cardinal = contains({results_fprime_choice_all(:).sessionStr},'Gr') & ...
                strcmp({results_fprime_choice_all(:).task},'cardinal') & ...
                cell2mat({results_fprime_choice_all(:).cohr_level}) == 15;

idx_gremlin_oblique = contains({results_fprime_choice_all(:).sessionStr},'Gr') & ...
                strcmp({results_fprime_choice_all(:).task},'oblique') & ...
                cell2mat({results_fprime_choice_all(:).cohr_level}) == 15;



tmp = {results_fprime_choice_all(idx_rolo_cardinal).dprime_stimulus};
dprime_stimulus_cardinal_rolo = cat(1,tmp{:});

tmp = {results_fprime_choice_all(idx_rolo_oblique).dprime_stimulus};
dprime_stimulus_oblique_rolo = cat(1,tmp{:});

tmp = {results_fprime_choice_all(idx_gremlin_cardinal).dprime_stimulus};
dprime_stimulus_cardinal_gremlin = cat(1,tmp{:});

tmp = {results_fprime_choice_all(idx_gremlin_oblique).dprime_stimulus};
dprime_stimulus_oblique_gremlin = cat(1,tmp{:});


tmp = {results_fprime_choice_all(idx_rolo_cardinal).dprime_stimulus_sign};
dprime_stimulus_sign_cardinal_rolo = cat(1,tmp{:});

tmp = {results_fprime_choice_all(idx_rolo_oblique).dprime_stimulus_sign};
dprime_stimulus_sign_oblique_rolo = cat(1,tmp{:});

tmp = {results_fprime_choice_all(idx_gremlin_cardinal).dprime_stimulus_sign};
dprime_stimulus_sign_cardinal_gremlin = cat(1,tmp{:});

tmp = {results_fprime_choice_all(idx_gremlin_oblique).dprime_stimulus_sign};
dprime_stimulus_sign_oblique_gremlin = cat(1,tmp{:});

%%% compute index of stimulus tuning bias
idx_stimulus_bias_rolo = (median(dprime_stimulus_cardinal_rolo,'omitnan') - median(dprime_stimulus_oblique_rolo,'omitnan')) / ...
                            (median(dprime_stimulus_cardinal_rolo,'omitnan') + median(dprime_stimulus_oblique_rolo,'omitnan'));
idx_signbias_cardinal_rolo =  median(dprime_stimulus_sign_cardinal_rolo,'omitnan') / std(dprime_stimulus_sign_cardinal_rolo,'omitnan');
idx_signbias_oblique_rolo  =  median(dprime_stimulus_sign_oblique_rolo,'omitnan') / std(dprime_stimulus_sign_oblique_rolo,'omitnan');


idx_stimulus_bias_gremlin = (median(dprime_stimulus_cardinal_gremlin,'omitnan') - median(dprime_stimulus_oblique_gremlin,'omitnan')) / ...
                            (median(dprime_stimulus_cardinal_gremlin,'omitnan') + median(dprime_stimulus_oblique_gremlin,'omitnan'));
idx_signbias_cardinal_gremlin =  median(dprime_stimulus_sign_cardinal_gremlin,'omitnan') / std(dprime_stimulus_sign_cardinal_gremlin,'omitnan');
idx_signbias_oblique_gremlin  =  median(dprime_stimulus_sign_oblique_gremlin,'omitnan') / std(dprime_stimulus_sign_oblique_gremlin,'omitnan');

%% method 4: sample synthetic neuron with vonmises distribution until the global statistic is matched. no need to match shape of disribution

to_sample_bPF   = 0;
idx             = cell2mat({results_dprime(:).b_PF}) == to_sample_bPF; 

theta           = results_dprime(idx).projective_fields.phi_x;

dprime_stimulus_cardinal_model          = results_dprime(idx).dprime_cardinal;
dprime_stimulus_oblique_model           = results_dprime(idx).dprime_oblique;
dprime_stimulus_sign_cardinal_model     = results_dprime(idx).dprime_cardinal_sign;
dprime_stimulus_sign_oblique_model      = results_dprime(idx).dprime_oblique_sign;

b_sample_list = [-3:0.1:3];
tao_cardinal_list = [-3:0.1:3];
tao_oblique_list  = [-3:0.1:3];
nSample = 64;
K_bPF = numel(b_sample_list);
K_tao = numel(tao_cardinal_list);
% b_PF_list = [-1:0.4:0.8];
% tao_list = [-0.1:0.05:0.1];
[idx_stimulus_bias, idx_signbias_cardinal, idx_signbias_oblique] = deal(zeros(K_bPF,K_tao,K_tao));
[b_sample_all, tao_cardinal_all, tao_oblique_all] = deal(zeros(K_bPF,K_tao,K_tao));
for i = 1:numel(b_sample_list)
    for j1 = 1:numel(tao_cardinal_list)
        for j2 = 1:numel(tao_oblique_list)
        b_sample = b_sample_list(i);
       % tao = tao_list(j);
        tao_cardinal  = tao_cardinal_list(j1);
        tao_oblique   = tao_oblique_list(j2); 
            % b_PF = 1.8;
            % tao  = 0.1;
            [idx_sample] = util_it.sampling_vonMises(theta, b_sample, tao_cardinal, tao_oblique, nSample);
            
            dprime_cardinal_sampled = dprime_stimulus_cardinal_model(idx_sample);
            
            dprime_oblique_sampled  = dprime_stimulus_oblique_model(idx_sample);
            
            dprime_cardnal_sign_sampled = dprime_stimulus_sign_cardinal_model(idx_sample);
            
            dprime_oblique_sign_sampled = dprime_stimulus_sign_oblique_model(idx_sample);
            
            idx_stimulus_bias(i,j1, j2) = (median(dprime_cardinal_sampled) - median(dprime_oblique_sampled)) / ...
                                    (median(dprime_cardinal_sampled) + median(dprime_oblique_sampled));
            idx_signbias_cardinal(i,j1, j2) = median(dprime_cardnal_sign_sampled) / std(dprime_cardnal_sign_sampled);
            idx_signbias_oblique(i,j1, j2)  = median(dprime_oblique_sign_sampled) / std(dprime_oblique_sign_sampled);
            
            b_sample_all(i,j1,j2)            = b_sample;
            tao_cardinal_all(i,j1,j2)    = tao_cardinal;
            tao_oblique_all(i,j1,j2)     = tao_oblique; 
        end
    end
end

%%% find the combination of sampling parameter that makes the global statistic closet to empirical data
params_all          = [b_sample_all(:), tao_cardinal_all(:), tao_oblique_all(:)];
dat_model_all       = [idx_stimulus_bias(:), idx_signbias_cardinal(:), idx_signbias_oblique(:)];
dat_rolo_all        = [idx_stimulus_bias_rolo, idx_signbias_cardinal_rolo, idx_signbias_oblique_rolo];
dat_gremlin_all     = [idx_stimulus_bias_gremlin, idx_signbias_cardinal_gremlin, idx_signbias_oblique_gremlin];

idx_rolo            = knnsearch(dat_model_all, dat_rolo_all, 'K', 1);
idx_gremlin         = knnsearch(dat_model_all, dat_gremlin_all, 'K', 1);

%dat_rolo_all
matched_dat_rolo    = dat_model_all(idx_rolo,:);
params_rolo         = params_all(idx_rolo,:);

b_sample_rolo       = params_rolo(1); tao_cardinal_rolo = params_rolo(2); tao_oblique_rolo = params_rolo(3);
idx_sample_rolo     = util_it.sampling_vonMises(theta, b_sample_rolo, tao_cardinal_rolo, tao_oblique_rolo, nSample);

%dat_gremlin_all
matched_dat_gremlin = dat_model_all(idx_gremlin,:);
params_gremlin  = params_all(idx_gremlin,:);

b_sample_gremlin    = params_gremlin(1); tao_cardinal_gremlin = params_gremlin(2); tao_oblique_gremlin = params_gremlin(3);
idx_sample_gremlin  = util_it.sampling_vonMises(theta, b_sample_gremlin, tao_cardinal_gremlin, tao_oblique_gremlin, nSample);



results_rolo.idx_stimulus_bias              = idx_stimulus_bias_rolo;
results_rolo.idx_signbias_cardinal          = idx_signbias_cardinal_rolo;
results_rolo.idx_signbias_oblique           = idx_signbias_oblique_rolo;
results_rolo.idx_stimulus_bias_matched      = matched_dat_rolo(1);
results_rolo.idx_signbias_cardinal_matched  = matched_dat_rolo(2);
results_rolo.idx_signbias_oblique_matched   = matched_dat_rolo(3);
results_rolo.idx_sample                     = idx_sample_rolo;
results_rolo.b_sample                       = b_sample_rolo;
results_rolo.tao_cardinal                   = tao_cardinal_rolo;
results_rolo.tao_oblique                    = tao_oblique_rolo;

%%% save idx of sampled synthetic neuron for future use
save_folder = '../../results/filtered_neuron_synthetic';
save_name = fullfile(save_folder,sprintf('sampled_subset_empirical_rolo_nTotal_%d_nSample_%d_b_PF_%s',...
                    nNeuron, nSample, strrep(sprintf('%.2f', to_sample_bPF), '.', '_')));
save(save_name,'results_rolo')



results_gremlin.idx_stimulus_bias              = idx_stimulus_bias_gremlin;
results_gremlin.idx_signbias_cardinal          = idx_signbias_cardinal_gremlin;
results_gremlin.idx_signbias_oblique           = idx_signbias_oblique_gremlin;
results_gremlin.idx_stimulus_bias_matched      = matched_dat_gremlin(1);
results_gremlin.idx_signbias_cardinal_matched  = matched_dat_gremlin(2);
results_gremlin.idx_signbias_oblique_matched   = matched_dat_gremlin(3);
results_gremlin.idx_sample                     = idx_sample_gremlin;
results_gremlin.b_sample                       = b_sample_gremlin;
results_gremlin.tao_cardinal                   = tao_cardinal_gremlin;
results_gremlin.tao_oblique                    = tao_oblique_gremlin;
save_name = fullfile(save_folder,sprintf('sampled_subset_empirical_gremlin_nTotal_%d_nSample_%d_b_PF_%s',...
                    nNeuron, nSample, strrep(sprintf('%.2f', to_sample_bPF), '.', '_')));
save(save_name,'results_gremlin')

%% Unsucessful methods
%%% Method 1: sample based on dprime_cardinal and dprime_oblique separately 
doThis = 0;
if doThis
    sample_options.doNorm       = true;
    
    sample_options.doPlot       = true;
    sample_options.nSample      = 128;
    
    
    if sample_options.doNorm
        data_model_cardinal     = dprime_stimulus_sign_cardinal_model_norm';
        data_model_oblique      = dprime_stimulus_sign_oblique_model_norm';
        data_rolo_cardinal      = dprime_stimulus_sign_cardinal_rolo_norm;
        data_rolo_oblique       = dprime_stimulus_sign_oblique_rolo_norm;
        data_gremlin_cardinal   = dprime_stimulus_sign_cardinal_gremlin_norm;
        data_gremlin_oblique    = dprime_stimulus_sign_oblique_gremlin_norm;
    else
        data_model_cardinal     = dprime_stimulus_sign_cardinal_model';
        data_model_oblique      = dprime_stimulus_sign_oblique_model';
        data_rolo_cardinal      = dprime_stimulus_sign_cardinal_rolo;
        data_rolo_oblique       = dprime_stimulus_sign_oblique_rolo;
        data_gremlin_cardinal   = dprime_stimulus_sign_cardinal_gremlin;
        data_gremlin_oblique    = dprime_stimulus_sign_oblique_gremlin;
    end
    [~,idx_sample_cardinal_rolo] = sample_match_distribution(data_model_cardinal, data_rolo_cardinal, sample_options);
    [~,idx_sample_oblique_rolo]  = sample_match_distribution(data_model_oblique, data_rolo_oblique, sample_options);
    
    [~,idx_sample_cardinal_gremlin] = sample_match_distribution(data_model_cardinal, data_gremlin_cardinal, sample_options);
    [~,idx_sample_oblique_gremlin]  = sample_match_distribution(data_model_oblique, data_gremlin_oblique, sample_options);
    
    idx_sample_rolo_interset = intersect(idx_sample_cardinal_rolo,idx_sample_oblique_rolo);
    idx_sample_rolo_union     = union(idx_sample_cardinal_rolo,idx_sample_oblique_rolo);
    
    idx_sample_gremlin_interset = intersect(idx_sample_cardinal_gremlin,idx_sample_oblique_gremlin);
    idx_sample_gremlin_union    = union(idx_sample_cardinal_gremlin,idx_sample_oblique_gremlin);
    
    
    
    
    figure;
    subplot(2,2,1); hold on
    cdfplot(data_model_cardinal);cdfplot(data_rolo_cardinal); 
    cdfplot(data_model_cardinal(idx_sample_rolo_interset));
    cdfplot(data_model_cardinal(idx_sample_rolo_union));
    subplot(2,2,2); hold on
    histogram(data_model_cardinal,'Normalization','probability');
    histogram(data_rolo_cardinal,'Normalization','probability'); 
    histogram(data_model_cardinal(idx_sample_rolo_interset), 'Normalization','probability');
    histogram(data_model_cardinal(idx_sample_rolo_union), 'Normalization','probability');
    
    subplot(2,2,3); hold on
    cdfplot(data_model_oblique);cdfplot(data_rolo_oblique); 
    cdfplot(data_model_oblique(idx_sample_rolo_interset));
    cdfplot(data_model_oblique(idx_sample_rolo_union));
    subplot(2,2,4); hold on
    histogram(data_model_oblique,'Normalization','probability');
    histogram(data_rolo_oblique,'Normalization','probability'); 
    histogram(data_model_oblique(idx_sample_rolo_interset), 'Normalization','probability');
    histogram(data_model_oblique(idx_sample_rolo_union), 'Normalization','probability');
    
    
    figure;
    subplot(2,2,1); hold on
    cdfplot(data_model_cardinal);cdfplot(data_gremlin_cardinal); 
    cdfplot(data_model_cardinal(idx_sample_gremlin_interset));
    cdfplot(data_model_cardinal(idx_sample_gremlin_union));
    subplot(2,2,2); hold on
    histogram(data_model_cardinal,'Normalization','probability');
    histogram(data_rolo_cardinal,'Normalization','probability'); 
    histogram(data_model_cardinal(idx_sample_gremlin_interset), 'Normalization','probability');
    histogram(data_model_cardinal(idx_sample_gremlin_union), 'Normalization','probability');
    
    subplot(2,2,3); hold on
    cdfplot(data_model_oblique);cdfplot(data_gremlin_oblique); 
    cdfplot(data_model_oblique(idx_sample_gremlin_interset));
    cdfplot(data_model_oblique(idx_sample_gremlin_union));
    subplot(2,2,4); hold on
    histogram(data_model_oblique,'Normalization','probability');
    histogram(data_rolo_oblique,'Normalization','probability'); 
    histogram(data_model_oblique(idx_sample_gremlin_interset), 'Normalization','probability');
    histogram(data_model_oblique(idx_sample_gremlin_union), 'Normalization','probability');
end
%%% Method 2: matching the difference between two distribution

doThis = 0;
if doThis
    sample_options.doNorm       = true;
    
    sample_options.doPlot       = true;
    sample_options.nSample      = 128;
    
    
    if sample_options.doNorm
        data_model      = dprime_stimulus_sign_diff_model;
        data_rolo       = dprime_stimulus_sign_diff_rolo;
        data_gremlin    = dprime_stimulus_sign_diff_gremlin;
    else
        data_model      = dprime_stimulus_sign_diff_model_norm;
        data_rolo       = dprime_stimulus_sign_diff_rolo_norm;
        data_gremlin    = dprime_stimulus_sign_diff_gremlin_norm;
    end
    [~,idx_sample_rolo]         = sample_match_distribution(data_model, data_rolo, sample_options);
    [~,idx_sample_gremlin]      = sample_match_distribution(data_model, data_gremlin, sample_options);
    
    
    
    
    figure;
    subplot(2,2,1); hold on
    cdfplot(data_model_cardinal);cdfplot(data_rolo_cardinal); 
    cdfplot(data_model_cardinal(idx_sample_rolo));
    cdfplot(data_model_cardinal(idx_sample_rolo));
    subplot(2,2,2); hold on
    histogram(data_model_cardinal,'Normalization','probability');
    histogram(data_rolo_cardinal,'Normalization','probability'); 
    histogram(data_model_cardinal(idx_sample_rolo), 'Normalization','probability');
    histogram(data_model_cardinal(idx_sample_rolo), 'Normalization','probability');
    
    subplot(2,2,3); hold on
    cdfplot(data_model_oblique);cdfplot(data_rolo_oblique); 
    cdfplot(data_model_oblique(idx_sample_gremlin));
    cdfplot(data_model_oblique(idx_sample_gremlin));
    subplot(2,2,4); hold on
    histogram(data_model_oblique,'Normalization','probability');
    histogram(data_rolo_oblique,'Normalization','probability'); 
    histogram(data_model_oblique(idx_sample_gremlin), 'Normalization','probability');
    histogram(data_model_oblique(idx_sample_gremlin), 'Normalization','probability');
end
%%% method 3: sample two dimensional dprime value joitly
doThis = 0;
if doThis
    sample_options.doNorm       = true;
    
    sample_options.nSample      = 128;
    
    
    if sample_options.doNorm
        data_model_cardinal     = dprime_stimulus_sign_cardinal_model_norm';
        data_model_oblique      = dprime_stimulus_sign_oblique_model_norm';
        data_rolo_cardinal      = dprime_stimulus_sign_cardinal_rolo_norm;
        data_rolo_oblique       = dprime_stimulus_sign_oblique_rolo_norm;
        data_gremlin_cardinal   = dprime_stimulus_sign_cardinal_gremlin_norm;
        data_gremlin_oblique    = dprime_stimulus_sign_oblique_gremlin_norm;
    else
        data_model_cardinal     = dprime_stimulus_sign_cardinal_model';
        data_model_oblique      = dprime_stimulus_sign_oblique_model';
        data_rolo_cardinal      = dprime_stimulus_sign_cardinal_rolo;
        data_rolo_oblique       = dprime_stimulus_sign_oblique_rolo;
        data_gremlin_cardinal   = dprime_stimulus_sign_cardinal_gremlin;
        data_gremlin_oblique    = dprime_stimulus_sign_oblique_gremlin;
    end
    
    
    data_model      = [data_model_cardinal, data_model_oblique];
    data_rolo       = [data_rolo_cardinal, data_rolo_oblique];
    data_gremlin    = [data_gremlin_cardinal, data_gremlin_oblique]; 
    [~, idx_sample_rolo] = sample_match_distribution_knn(data_model, data_rolo, sample_options);
    [~, idx_sample_gremlin] = sample_match_distribution_knn(data_model, data_gremlin, sample_options);
    
    figure
    subplot(2,2,1); hold on
    cdfplot(data_model_cardinal);cdfplot(data_rolo_cardinal); 
    cdfplot(data_model_cardinal(idx_sample_rolo));
    legend('Synthetic','Empirical','Sampled'); box off
    title('Rolo cardinal')
    set(gca,'fontsize',18); xlabel('dprime')

    subplot(2,2,2); hold on
    cdfplot(data_model_oblique);cdfplot(data_rolo_oblique); 
    cdfplot(data_model_oblique(idx_sample_rolo));
    legend('Synthetic','Empirical','Sampled'); box off
    title('Rolo oblique')
    set(gca,'fontsize',18); xlabel('dprime')

    subplot(2,2,3); hold on
    cdfplot(data_model_cardinal);cdfplot(data_gremlin_cardinal); 
    cdfplot(data_model_cardinal(idx_sample_gremlin));
    legend('Synthetic','Empirical','Sampled'); box off
    title('Gremlin cardinal')
    set(gca,'fontsize',18); xlabel('dprime')

    subplot(2,2,4); hold on
    cdfplot(data_model_oblique);cdfplot(data_gremlin_oblique); 
    cdfplot(data_model_oblique(idx_sample_gremlin));
    legend('Synthetic','Empirical','Sampled'); box off
    title('Gremlin oblique')
    set(gca,'fontsize',18); xlabel('dprime')
end
%% helper functions


function phi = sample_phi_signed(tao, b_PF, nNeuron)
%%%%% This code was adapted from C_projection. Other than sampling
%%%%% orientation with a cardnal or oblique bias (controlled by b_PF), it
%%%%% also considers a sign bias (controlled by tao) which controls whether
%%%%% more neurons prefer positive or negative orientation

theta = [0:0.01:180] / 180 * pi;
%%%% 
if b_PF > 0
    %%% overrepresent cardinal orientation
    y =  (exp(b_PF * cos(2 * theta))+exp(b_PF * cos(2*(theta - pi/2))));
else
    %%%% overrepresent oblique orientation
    y =  (exp(b_PF * cos(2 * (theta - pi/4)))+exp(b_PF * cos(2*(theta - 3 *pi/4))));
end
%%%%% add sign bias
y = y .* exp(-theta * tao);

probabilities = y / sum(y); % Normalize to sum to 1

% Create cumulative distribution function (CDF)
cumdf =  cumsum(probabilities); 
cumdf(1) = 0; % Add 0 at the start

% Generate uniform random samples
linearProb = linspace(0,1,nNeuron);
phi = interp1(cumdf,theta, linearProb, 'linear','extrap');


end

function [sampled_values, idx_sample] = sample_match_distribution_knn(list_original, list_target, sample_options)


nTarget = size(list_target, 1);
nOriginal = size(list_original, 1);
nSample  = sample_options.nSample;
% if nTarget > nOriginal
%     error('Target list has more points than original list. Unique matching is not possible.');
% end

% % Build KD-tree for efficient NN search
% Mdl = createns(list_origin
idx_target      = randsample(nTarget,nSample);
target_values   = list_target(idx_target,:);

idx_sample      = knnsearch(list_original, target_values);
idx_sample      = unique(idx_sample);
sampled_values  = list_original(idx_sample, :);

% % Initialize
% used_idx = false(nOriginal, 1);
% idx_sample = zeros(nSample, 1);
% nK  = 5;
% for i = 1:nSample
%     % Find the nearest available neighbor
%     [idx_all, ~] = knnsearch(list_original, target_values(i,:), 'K', nK);
% 
%     % Select the first unused neighbor
%     for j = 1:nK
%         if ~used_idx(idx_all(j))
%             idx_sample(i) = idx_all(j);
%             used_idx(idx_all(j)) = true;
%             break;
%         end
%     end
% end
% idx_sample      = unique(nonzeros(idx_sample));
% sampled_values = list_original(idx_sample, :);
end


function [sampled_values, idx_sample_real] = sample_match_distribution(list_original, list_target, options)


nOriginal       = numel(list_original); 
sorted_original = sort(list_original);
%p_original      = linspace(0, 1, nOriginal);

% Step 1: Get target quantiles
quantile_edges = linspace(0, 1, options.nSample);
target_quantiles = quantile(list_target, quantile_edges);

% % Step 2: Inverse-CDF mapping from quantiles to list_original values
 p_original = linspace(0, 1, nOriginal);
% values_from_original = interp1(p_original, sorted_original, target_quantiles, 'linear', 'extrap');

idx_sample = round(nOriginal * interp1(sorted_original, p_original, target_quantiles, 'linear', 'extrap'));
idx_sample(idx_sample <= 0) = 1;
idx_sample(idx_sample > nOriginal) = nOriginal;
idx_sample = unique(idx_sample);
sampled_values = sorted_original(idx_sample);

[~, idx_sample_real] = ismember(sampled_values, list_original);
    if options.doPlot
        figure;
        subplot(1,2,1)
        cdfplot(list_original); hold on
        cdfplot(list_target); 
        cdfplot(sampled_values)
        subplot(1,2,2)
        histogram(list_original,'Normalization','probability');
        hold on; histogram(list_target,'Normalization','probability'); 
        histogram(sampled_values,'Normalization','probability')
    end
end
% %%
% % Sample input

% 
% % Step 1: Estimate empirical CDF of listB
% cdf_B_x = sort(listB);
% cdf_B_p = linspace(0, 1, numel(cdf_B_x));
% 
% % Step 2: Get quantile targets from listB
% target_probs = rand(numel(listB), 1);  % Uniform samples
% target_values = interp1(cdf_B_p, cdf_B_x, target_probs, 'linear', 'extrap');
% 
% % Step 3: Build inverse CDF (quantile function) of listA
% sorted_A = sort(listA);
% p_A = linspace(0, 1, numel(sorted_A));
% 
% % Step 4: Find listA values that correspond to listB's quantiles
% sampled_from_A = interp1(p_A, sorted_A, ...
%                          interp1(cdf_B_x, cdf_B_p, target_values, 'linear', 'extrap'), ...
%                          'linear', 'extrap');
% 
% % `sampled_from_A` is now a set of values from listA that approximate listB's distribution
