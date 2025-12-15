clear all
clc
close all
%%%%%% This script checks the number of zero-coherence trials and number of
%%%%%% neurons for each session, to make sure we still can pass the
%%%%%% criterion after separating zero-coherence trials two halves by eye
%%%%%% position
%%
global bpGlobal
bpGratingFCGlobal;
session_list_ro = bpGlobal.rolo.session_list.switching;
session_list_gr = bpGlobal.gremlin.session_list.interleaved_real;

session_list_all = [session_list_ro;session_list_gr];
%% check the trial number of zero-coherence trials and neurons
neuron_filter_name = 'coef1_hVis2_FR1';
load(sprintf('/Users/liushizhao/projects_local/bpGratingEx/results/filter_neuron/%s/filtered_neurons_%s',neuron_filter_name, neuron_filter_name));

results_numbers = struct();
for n = 1:numel(session_list_all)
    sessionStr      = session_list_all{n};
    subjectCode     = sessionStr(1:2);

    data_folder = sprintf('/Users/liushizhao/projectData_local/probinf_data/extractedData/%s/bpGratingFC', subjectCode);
    dlist = dir(fullfile(data_folder, [sessionStr,'*.mat'] ));
    load(fullfile(data_folder, dlist(1).name),'stimulus','behavior');

    nZero_cardinal  = sum(strcmp(stimulus.taskType,'C') & stimulus.signal' == 0 & behavior.completeOnline');
    nZero_oblique   = sum(strcmp(stimulus.taskType,'O') & stimulus.signal' == 0 & behavior.completeOnline');

    %%%% number of neurons
    idx_session = strcmp({filtered_neurons(:).sessionStr},sessionStr);
    nNeuron = filtered_neurons(idx_session).nNeuron_kept;

    results_numbers(n).sessionStr = sessionStr;
    results_numbers(n).nNeuron = nNeuron;
    results_numbers(n).nZero_cardinal = nZero_cardinal;
    results_numbers(n).nZero_cardinal_half = floor(nZero_cardinal/2);
    results_numbers(n).nZero_oblique = nZero_oblique;
    results_numbers(n).nZero_oblique_half = floor(nZero_oblique/2);

end
%% see how many sessions did not pass the criterion
%%% T_zero_half > N + 4
figure;
subplot(2,2,1)
idx = ismember(session_list_all, session_list_ro);
nNeuron_all = cell2mat({results_numbers(idx).nNeuron});
nZero_cardinal_half = cell2mat({results_numbers(idx).nZero_cardinal_half});
nZero_oblique_half  = cell2mat({results_numbers(idx).nZero_oblique_half});
nZero_half = min([nZero_cardinal_half; nZero_oblique_half],[],1);
histogram(nZero_half - nNeuron_all - 4);
subplot(2,2,3)
plot(nZero_half - 4, '-o'); ylim([0,50])
yyaxis right
plot(nNeuron_all, '-o');ylim([0,50])

subplot(2,2,2)
idx = ismember(session_list_all, session_list_gr);
nNeuron_all = cell2mat({results_numbers(idx).nNeuron});
nZero_cardinal_half = cell2mat({results_numbers(idx).nZero_cardinal_half});
nZero_oblique_half  = cell2mat({results_numbers(idx).nZero_oblique_half});
nZero_half = min([nZero_cardinal_half; nZero_oblique_half],[],1);
histogram(nZero_half - nNeuron_all - 4);
subplot(2,2,4)
plot(nZero_half - 4, '-o');ylim([0,max(nNeuron_all) + 5])
yyaxis right
plot(nNeuron_all, '-o');ylim([0,max(nNeuron_all) + 5 ])


%% number of neuron that meets threshold for each session
figure
subplot(1,2,1)
idx = ismember(session_list_all, session_list_ro);
nZero_cardinal_half = cell2mat({results_numbers(idx).nZero_cardinal_half});
nZero_oblique_half  = cell2mat({results_numbers(idx).nZero_oblique_half});
nZero_half = min([nZero_cardinal_half; nZero_oblique_half],[],1);
plot(nZero_half - 4, '-o')
yyaxis


subplot(1,2,2)
idx = ismember(session_list_all, session_list_gr);
nZero_cardinal_half = cell2mat({results_numbers(idx).nZero_cardinal_half});
nZero_oblique_half  = cell2mat({results_numbers(idx).nZero_oblique_half});
nZero_half = min([nZero_cardinal_half; nZero_oblique_half],[],1);
plot(nZero_half - 4, '-o')