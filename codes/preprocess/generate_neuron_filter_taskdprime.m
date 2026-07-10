clear all
clc
close all
%%%%%% This script generate neuron filter based on stimulus tuning of
%%%%%% neurons. We separate the population into two groups of equal size
%%%%%% based on their stimulus dprime
%% prepare data
global bpGlobal
bpGratingFCGlobal;

%%% data path and session list
%%%%%%%%%%%%%%%%%%%%%%%%%% specify session list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% session list of rolo
session_list_ro = bpGlobal.rolo.session_list;
session_list_ro_all = session_list_ro.switching;

%%%% session list of gremlin
session_list_gr = bpGlobal.gremlin.session_list;
session_list_gr_all = session_list_gr.interleaved_real;

% session_rolo_cardinal       = session_list_ro.cardinal;
% session_rolo_oblique        = session_list_ro.oblique;
% session_gremlin_cardinal    = session_list_gr.cardinal_noError;
% session_gremlin_oblique     = session_list_gr.oblique_use;

session_list_all    = [session_list_ro_all; session_list_gr_all];
nSession            = numel(session_list_all);

% cardinal_sessions           = [session_rolo_cardinal;session_gremlin_cardinal];
% oblique_sessions            = [session_rolo_oblique;session_gremlin_oblique];

%% generate neuron filter
%%% highC_highO, lowC_highO
%%% highC_lowO, lowC_lowO
%%%% cross_combination: highC_lowO_lowC_highO
%%%%% same_combination: highC_highO_lowC_lowO
doThis = 0;
if doThis
    dprime_option   = 'highC_highO_lowC_lowO';
    session_type    = 'passiveViewing';
    filter_name     = sprintf('%s_%s',dprime_option,session_type); 
    criterion_message = [''];
    switch session_type
        case 'mainTask'
            load('../../results/neural/stimulus_dprime_twoanimals_interleaved.mat');
            % switch dprime_option
            %     case 'high'
            %         criterion_message   = ['Keep half of the neurons of high stimulus dprime.\n'...
            %             sprintf('dprime was computed with this version: %s',results_fprime_choice_all(1).versionName) ];
            %         filter_name         = 'higher_stimulusdprime';
            %     case 'low'
            %         criterion_message   = ['Keep half of the neurons of low stimulus dprime.\n'...
            %             sprintf('dprime was computed with this version: %s',results_fprime_choice_all(1).versionName) ];
            %         filter_name         = 'lower_stimulusdprime';
            % end
        case 'passiveViewing'
            load('../../results/neural/stimulus_dprime_twoanimals_interleaved_passive.mat');
            % switch dprime_option
            %     case 'high'
            %         criterion_message   = ['Keep half of the neurons of high stimulus dprime.\n'...
            %             'dprime was computed using passive viewing session.\n'...
            %             sprintf('dprime was computed with this version: %s',results_fprime_choice_passive_all(1).versionName) ];
            %         filter_name         = 'higher_stimulusdprime_passive';
            %     case 'low'
            %         criterion_message   = ['Keep half of the neurons of low stimulus dprime.\n'...
            %             'dprime was computed using passive viewing session.\n'...
            %             sprintf('dprime was computed with this version: %s',results_fprime_choice_passive_all(1).versionName) ];
            %         filter_name         = 'lower_stimulusdprime_passive';
            % end
    end

    savePath            = fullfile('../../results/filter_neuron',filter_name);
    txtName             = fullfile(savePath,sprintf('criterion_%s.txt',filter_name));
    saveName            = fullfile(savePath,sprintf('filtered_neurons_%s',filter_name));

    %%%%%%%%%%%%%%%% save criterion message in a .txt file   %%%%%%%%%%%%
    mkdir(savePath);
    fid = fopen(txtName,'wt');
    fprintf(fid, criterion_message);
    fclose(fid);

    %%%%%%%%%%%%%%%%%%%%%%%%% also keep some crierion in a struct  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    keep_options                    = struct();
    keep_options.unitOptions        = struct();
    keep_options.unitFilter_name    = filter_name;

    for i = 1:nSession

        %i = 1;
        sessionStr  = session_list_all{i};
        subjectCode = sessionStr(1:2);
        dateStr     = sessionStr(3:end);
    
        %%%%% get the dprime results
        switch session_type
            case 'mainTask'
                idx_base = strcmp({results_fprime_choice_all(:).sessionStr},sessionStr) & ...
                    cell2mat({results_fprime_choice_all(:).cohr_level}) > 0 & ...
                     cell2mat({results_fprime_choice_all(:).cohr_level}) <= 20;
                
                idx_cardinal    = idx_base & strcmp({results_fprime_choice_all(:).task}, 'cardinal');
                idx_oblique     = idx_base & strcmp({results_fprime_choice_all(:).task}, 'oblique');
                
                results_session_cardinal    = results_fprime_choice_all(idx_cardinal);
                results_session_oblique     = results_fprime_choice_all(idx_oblique);

            case 'passiveViewing'
      
                idx_base = strcmp({results_fprime_choice_passive_all(:).sessionStr},sessionStr) & ...
                    cell2mat({results_fprime_choice_passive_all(:).cohr_level}) > 0 & ...
                     cell2mat({results_fprime_choice_passive_all(:).cohr_level}) <= 20;

                idx_cardinal    = idx_base & strcmp({results_fprime_choice_passive_all(:).task}, 'cardinal');
                idx_oblique     = idx_base & strcmp({results_fprime_choice_passive_all(:).task}, 'oblique');
                
                results_session_cardinal    = results_fprime_choice_passive_all(idx_cardinal);
                results_session_oblique     = results_fprime_choice_passive_all(idx_oblique);
                
        end
        if isempty(results_session_cardinal) | isempty(results_session_oblique)

            filtered_neurons(i).sessionStr          = [];
            continue
        end

        fisher_norm_avg_cardinal = get_fisher_norm_avg(results_session_cardinal);
        fisher_norm_avg_oblique = get_fisher_norm_avg(results_session_oblique);
        %%%% first square the dprime to get fisher information, which should be a
        %%%% linear function of coherence level, then normalize by coherence level.
    
        idx_high_cardinal           = fisher_norm_avg_cardinal >= median(fisher_norm_avg_cardinal, 'omitnan');
        idx_low_cardinal            = fisher_norm_avg_cardinal <  median(fisher_norm_avg_cardinal, 'omitnan');

        idx_high_oblique            = fisher_norm_avg_oblique >= median(fisher_norm_avg_oblique, 'omitnan');
        idx_low_oblique             = fisher_norm_avg_oblique <  median(fisher_norm_avg_oblique, 'omitnan');


        switch dprime_option
            case 'highC_highO'
                idx_keep = idx_high_cardinal & idx_high_oblique;
            case 'lowC_highO'
                idx_keep = idx_low_cardinal & idx_high_oblique;
            case 'highC_lowO'
                idx_keep = idx_high_cardinal & idx_low_oblique;
            case 'lowC_lowO'
                idx_keep = idx_low_cardinal & idx_low_oblique;
            case 'highC_lowO_lowC_highO'
                idx_keep =  (idx_high_cardinal & idx_low_oblique) |  (idx_low_cardinal & idx_high_oblique);
            case 'highC_highO_lowC_lowO'
                idx_keep = (idx_high_cardinal & idx_high_oblique) | (idx_low_cardinal & idx_low_oblique);
        end
        if sum(idx_keep) == 0
            sum(idx_keep)
        end
        neuronIdx_kept.unitId       = results_session_cardinal(1).unitId(idx_keep);
        neuronIdx_kept.electrode    = results_session_cardinal(1).electrode(idx_keep);
        neuronIdx_kept.idOnChannel  = results_session_cardinal(1).idOnChannel(idx_keep);

        nNeuron_kept    = numel(neuronIdx_kept.unitId);
        %nNeuron_total   = numel(unitProperties);

        %%%%%%%%%%%%%%%%% save all information in a struct %%%%%%%%%%%%%%%%%
        filtered_neurons(i).sessionStr           = sessionStr;
        filtered_neurons(i).subjectCode          = subjectCode;
        filtered_neurons(i).dateStr              = dateStr;

        filtered_neurons(i).neuronIdx_kept       = neuronIdx_kept; % kept neurons
        filtered_neurons(i).nNeuron_kept         = nNeuron_kept;
        %filtered_neurons(i).nNeuron_total        = nNeuron_total;
        filtered_neurons(i).unitOptions          = keep_options;



    end

    validRows = ~cellfun(@isempty, {filtered_neurons(:).sessionStr});
    filtered_neurons = filtered_neurons(validRows);

    save(saveName,'filtered_neurons','keep_options');
end
%% subsample neuron for size control


doThis = 1;
if doThis
    neuron_filter_name_list =  {'highC_highO_passiveViewing';'highC_lowO_passiveViewing';'lowC_highO_passiveViewing';'lowC_lowO_passiveViewing';...
                    'highC_highO_mainTask';'highC_lowO_mainTask';'lowC_highO_mainTask';'lowC_lowO_mainTask';...
                    'highC_lowO_lowC_highO_passiveViewing';'highC_highO_lowC_lowO_passiveViewing';...
                    'highC_lowO_lowC_highO_mainTask';'highC_highO_lowC_lowO_mainTask'};
    for k = 1:numel(neuron_filter_name_list)
   % neuron_filter_name  = 'highC_highO_passiveViewing';
        neuron_filter_name = neuron_filter_name_list{k};
        savePath            = fullfile('../../results/filter_neuron',neuron_filter_name);
        loadName            = fullfile(savePath,sprintf('filtered_neurons_%s',neuron_filter_name));
        load(loadName)
        nBootstrap          = 1000;
        %%%%%%% these numbers are from generate_neuron_sizeControl, version = coef1_hVis2_FR1
        switch neuron_filter_name
            case {'highC_highO_passiveViewing';'highC_lowO_passiveViewing';'lowC_highO_passiveViewing';'lowC_lowO_passiveViewing';...
                    'highC_highO_mainTask';'highC_lowO_mainTask';'lowC_highO_mainTask';'lowC_lowO_mainTask'}
                nSub_rolo       = 4;
                nSub_gremlin    = 52/4;
            case {'highC_lowO_lowC_highO_passiveViewing';'highC_highO_lowC_lowO_passiveViewing';...
                    'highC_lowO_lowC_highO_mainTask';'highC_highO_lowC_lowO_mainTask'}
                nSub_rolo       = 7;
                nSub_gremlin    = 52/2;
        end
    
        filtered_neurons_new_rolo = get_neuron_sizeControl(filtered_neurons, nSub_rolo, session_list_ro_all, nBootstrap);
       
        filtered_neurons_new_gremlin = get_neuron_sizeControl(filtered_neurons, nSub_gremlin, session_list_gr_all, nBootstrap);
    
    
        % filtered_neurons_new_rolo_interleaved    = get_neuron_sizeControl(filtered_neurons, nSub_rolo_interleaved, session_rolo_interleaved, nBootstrap);
        % filtered_neurons_new_gremlin_interleaved = get_neuron_sizeControl(filtered_neurons, nSub_gremlin_interleaved, session_gremlin_interleaved, nBootstrap);
        %
        filtered_neurons_new = [filtered_neurons_new_rolo,filtered_neurons_new_gremlin];
    
    
    
        filtered_neurons = filtered_neurons_new;
    
        savepath = sprintf('../../results/filter_neuron/%s_sizeControl',neuron_filter_name);
        mkdir(savepath)
        saveName = fullfile(savepath, sprintf('filtered_neurons_%s_sizeControl', neuron_filter_name));
    
        save(saveName, 'filtered_neurons','keep_options');
    end
end


%%
function fisher_norm_avg = get_fisher_norm_avg(results_session)

dprime_stimulus = {results_session(:).dprime_stimulus};
cohr            = cell2mat({results_session(:).cohr_level});
fisher_norm     = arrayfun(@(k)dprime_stimulus{k} .^ 2 / cohr(k), [1:numel(cohr)], ...
    'UniformOutput',false);
fisher_norm_avg = mean(cat(2,fisher_norm{:}),2,'omitnan');
end