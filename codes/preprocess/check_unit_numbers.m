clear all
clc
close all
%%% Number of units, single units and multiunits
global bpGlobal
bpGratingFCGlobal()
%%
session_list_rolo = bpGlobal.rolo.session_list.switching;
session_list_gremlin = bpGlobal.gremlin.session_list.interleaved_real;

snr_thres = 2.5;
%% 1. Number of units after spikesorting
subjectCode = 'Gr';
switch subjectCode
    case 'Ro'
        data_path = '/Users/liushizhao/projectData_local/probinf_data/extractedData/Ro/unitProperties';
        session_list = session_list_rolo;
    case 'Gr'
        data_path = '/Users/liushizhao/projectData_local/probinf_data/extractedData/Gr/unitProperties';
        session_list = session_list_gremlin;
end
nSession = numel(session_list);
[nTotal, nSingle, nMulti] = deal(zeros(nSession, 1));
for n = 1:numel(session_list)
    sessionStr = session_list{n};
    load(fullfile(data_path, sprintf('%s_unitProperties',sessionStr)));
    waveformSNR = cell2mat({unitProperties(:).waveformSNR});
    nTotal(n) = numel(unitProperties);
    nSingle(n) = sum(waveformSNR >= snr_thres);
    nMulti(n) = sum(waveformSNR < snr_thres);
end

fprintf('Inital spikesorting: Monkey %s, nUnits: %.2f +- %.2f. \n', subjectCode, mean(nTotal), std(nTotal));
fprintf('Inital spikesorting: Monkey %s, nSingleUnits: %.2f +- %.2f. \n', subjectCode, mean(nSingle), std(nSingle));
fprintf('Inital spikesorting: Monkey %s, nMultiUnits: %.2f +- %.2f. \n', subjectCode, mean(nMulti), std(nMulti));
%% 2. Number of units after applying some criteria
subjectCode = 'Gr';
%%% filter options coef1_hVis2_FR1_hVisOri2_FROri2 or coef1_hVis2_FR1
filter_name = 'coef1_hVis2_FR1_hVisOri2_FROri2'; 

switch subjectCode
    case 'Ro'
        data_path = '/Users/liushizhao/projectData_local/probinf_data/extractedData/Ro/unitProperties';
        session_list = session_list_rolo;
    case 'Gr'
        data_path = '/Users/liushizhao/projectData_local/probinf_data/extractedData/Gr/unitProperties';
        session_list = session_list_gremlin;
end
load(sprintf('/Users/liushizhao/Documents/projects/bpGratingEx/results/filter_neuron/%s/filtered_neurons_%s',...
            filter_name, filter_name));

nSession = numel(session_list);
[nTotal, nSingle, nMulti] = deal(zeros(nSession, 1));
for n = 1:numel(session_list)
    sessionStr = session_list{n};
    load(fullfile(data_path, sprintf('%s_unitProperties',sessionStr)));


    chanLabel   = cell2mat({unitProperties(:).chanLabel});
    idOnChannel = cell2mat({unitProperties(:).idOnChannel});
    waveformSNR = cell2mat({unitProperties(:).waveformSNR});

    idx_session = strcmp({filtered_neurons(:).sessionStr}, sessionStr);

    neuronIdx_kept = filtered_neurons(idx_session).neuronIdx_kept;

    nKept = numel(neuronIdx_kept.electrode);
    waveform_kept = zeros(nKept,1);
    for k = 1:nKept
        idx = chanLabel == neuronIdx_kept.electrode(k) & idOnChannel == neuronIdx_kept.idOnChannel(k);
        waveform_kept(k) = waveformSNR(idx);
    end


    nTotal(n) = nKept;
    nSingle(n) = sum(waveform_kept >= snr_thres);
    nMulti(n)  = sum(waveform_kept < snr_thres);
end



fprintf('%s: Monkey %s, nUnits: %.2f +- %.2f. \n', filter_name, subjectCode, mean(nTotal), std(nTotal));
fprintf('%s: Monkey %s, nSingleUnits: %.2f +- %.2f. \n', filter_name, subjectCode, mean(nSingle), std(nSingle));
fprintf('%s: Monkey %s, nMultiUnits: %.2f +- %.2f. \n', filter_name, subjectCode, mean(nMulti), std(nMulti));