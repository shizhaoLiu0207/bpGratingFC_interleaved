clear all
clc
close all
%%%%% This script is used to generate subset of neurons with controled
%%%%% population size: one consistent size across all sessions in each
%%%%% epoch
%%
nBootstrap = 1000;
neuron_filter_name = 'coef1_hVis2_FR1';
size_option = 'sizeControl'; % sizeControl or half
save_filter_name   = [neuron_filter_name,'_interleaved_sizeControl'];

load(sprintf('/Users/liushizhao/Documents/projects/bpGratingEx/results/filter_neuron/%s/filtered_neurons_%s',neuron_filter_name, neuron_filter_name));
savepath = sprintf('../../results/filter_neuron/%s',save_filter_name);
mkdir(savepath)
saveName = fullfile(savepath, sprintf('filtered_neurons_%s', save_filter_name));

global bpGlobal
bpGratingFCGlobal;

%%%%% session list of rolo
session_list_ro = bpGlobal.rolo.session_list;
%session_list_ro_all = [session_list_ro.cardinal;session_list_ro.oblique;session_list_ro.switching];

%%%% session list of gremlin
session_list_gr = bpGlobal.gremlin.session_list;


% %%%  two interleaved epochs %%%%%
session_rolo_interleaved    = session_list_ro.switching;
session_gremlin_interleaved = session_list_gr.interleaved_real;

figure
subplot(1,2,1)
idx = ismember({filtered_neurons(:).sessionStr}, session_rolo_interleaved);
nUnit_kept = cell2mat({filtered_neurons(idx).nNeuron_kept});plot(nUnit_kept,'-o'); grid on; hold on
ylabel('Num. Neuron'); xlabel('Session'); set(gca,'fontsize', 20)
nUnit_sort = sort(nUnit_kept,'ascend');
plot(nUnit_sort,'-o')
title({'rolo interleaved';sprintf('min = %d, 2nd min = %d',nUnit_sort(1), nUnit_sort(2))})
subplot(1,2,2)
idx = ismember({filtered_neurons(:).sessionStr}, session_gremlin_interleaved);
nUnit_kept = cell2mat({filtered_neurons(idx).nNeuron_kept});plot(nUnit_kept,'-o'); grid on; hold on
ylabel('Num. Neuron'); xlabel('Session'); set(gca,'fontsize', 20)
nUnit_sort = sort(nUnit_kept,'ascend');
plot(nUnit_sort,'-o')
title({'gremlin interleaved';sprintf('min = %d, 2nd min = %d',nUnit_sort(1), nUnit_sort(2))})

switch neuron_filter_name
    case 'coef1_hVis2_FR1'

        nSub_rolo_interleaved = 14; %%%%% use the second minimum, the first minimum is 11
        nSub_gremlin_interleaved = 52;
    case 'coef1_hVis2_FR1_hVisOri2_FROri2'

        nSub_rolo_interleaved = 14; %%%%% just use the same as the other filter. but has to remove 6 sessions
        nSub_gremlin_interleaved = 38;
end

%%

filtered_neurons_new_rolo_interleaved    = get_neuron_sizeControl(filtered_neurons, nSub_rolo_interleaved, session_rolo_interleaved, nBootstrap);
filtered_neurons_new_gremlin_interleaved = get_neuron_sizeControl(filtered_neurons, nSub_gremlin_interleaved, session_gremlin_interleaved, nBootstrap);


filtered_neurons_new = [filtered_neurons_new_rolo_interleaved,filtered_neurons_new_gremlin_interleaved];

%%
filtered_neurons = filtered_neurons_new;
save(saveName, 'filtered_neurons','keep_options')
%% helper function 
function filtered_neurons_new = get_neuron_sizeControl(filtered_neurons, nSub, session_list, nBootstrap)
for n = 1:numel(session_list)
    
    idx = strcmp({filtered_neurons(:).sessionStr}, session_list{n});

   

    nKept = filtered_neurons(idx).nNeuron_kept;
    if nKept >= nSub
        if nchoosek(nKept,nSub) <= nBootstrap
            combination_all = nchoosek([1:nKept],nSub); 
        elseif nchoosek(nKept,nSub) <= 5 * nBootstrap
            combination_all = nchoosek([1:nKept],nSub); 
            combination_all = combination_all(1:nBootstrap,:); % each row is one combination
        else
            combination_all = zeros(nBootstrap,nSub);
            for k = 1:nBootstrap
                 combination_all(k,:) = randsample(nKept,nSub,'false');
            end
        end


        %%%%% Get the unitId, electrode, and idOnChannel of these combinations
        nCombine = size(combination_all,1);
        neuronIdx_kept_sub = struct();
        for k = 1:nCombine
            neuronIdx_kept_sub(k).unitId        = filtered_neurons(idx).neuronIdx_kept.unitId(combination_all(k,:));
            neuronIdx_kept_sub(k).electrode     = filtered_neurons(idx).neuronIdx_kept.electrode(combination_all(k,:));
            neuronIdx_kept_sub(k).idOnChannel   = filtered_neurons(idx).neuronIdx_kept.idOnChannel(combination_all(k,:));
        end
        filtered_neurons_new(n) = filtered_neurons(idx); % first copy everything
        filtered_neurons_new(n).nNeuron_kept = nSub;
        filtered_neurons_new(n).neuronIdx_kept = neuronIdx_kept_sub;
        
    end
    idx_remove_session = cellfun(@isempty,{filtered_neurons_new(:).sessionStr});

    filtered_neurons_new(idx_remove_session) = [];
   
end


end

