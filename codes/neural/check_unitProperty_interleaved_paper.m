clear all
clc
close all
%% 
filter_name = 'coef1_hVis2_FR1';
load(sprintf('../../results/filter_neuron/%s/filtered_neurons_%s',filter_name,filter_name));

global   bpGlobal ftsize
bpGratingFCGlobal();

session_list_gremlin    = bpGlobal.gremlin.session_list.interleaved_real;

session_list_rolo       = bpGlobal.rolo.session_list.switching;

%%
SNR_thres = 2.5;
subjectCode = 'Gr';

switch subjectCode
    case 'Ro'
        session_list = session_list_rolo;
    case 'Gr'
         session_list = session_list_gremlin;
end
unit_basic = struct();
for n = 1:numel(session_list)
    sessionStr = session_list{n};
    unitProperty_path = sprintf('/Users/liushizhao/projectData_local/probinf_data/extractedData/%s/unitProperties',sessionStr(1:2));
    load(fullfile(unitProperty_path,  sprintf('%s_unitProperties',sessionStr)));

    idx_session = strcmp({filtered_neurons(:).sessionStr}, sessionStr);
    neuronIdx_kept  = filtered_neurons(idx_session).neuronIdx_kept;
    nNeuron_kept    = filtered_neurons(idx_session).nNeuron_kept; 
    nNeuron_total   = filtered_neurons(idx_session).nNeuron_total;
    
    waveformSNR_neuronKept = zeros(1,nNeuron_kept);
    for i = 1:nNeuron_kept
        electrode   = neuronIdx_kept.electrode(i);
        idOnChannel      = neuronIdx_kept.idOnChannel(i);
        

        idx_neuron = cell2mat({unitProperties(:).chanLabel}) ==  electrode & ...
                     cell2mat({unitProperties(:).idOnChannel}) == idOnChannel;

        waveformSNR_neuronKept(i) = unitProperties(idx_neuron).waveformSNR;
    end
    waveformSNR_neuronAll = cell2mat({unitProperties(:).waveformSNR});


    unit_basic(n).sessionStr = sessionStr;
    unit_basic(n).nNeuron_total = nNeuron_total;
    unit_basic(n).nNeuron_kept  = nNeuron_kept;
    
    unit_basic(n).nSingleunit_kept = sum(waveformSNR_neuronKept >= SNR_thres);
    unit_basic(n).nMultiunit_kept  = sum(waveformSNR_neuronKept < SNR_thres);

    unit_basic(n).waveformSNR_neuronAll = waveformSNR_neuronAll;
    unit_basic(n).waveformSNR_neuronKept = waveformSNR_neuronKept;
end
waveformSNR_neuronAll = {unit_basic(:).waveformSNR_neuronAll};
waveformSNR_neuronAll = cat(2,waveformSNR_neuronAll{:});

waveformSNR_neuronKept = {unit_basic(:).waveformSNR_neuronKept};
waveformSNR_neuronKept = cat(2,waveformSNR_neuronKept{:});

nNeuron_total       = cell2mat({unit_basic(:).nNeuron_total}); 
nNeuron_kept       = cell2mat({unit_basic(:).nNeuron_kept}); 
nSingleunit_kept       = cell2mat({unit_basic(:).nSingleunit_kept}); 
nMultiunit_kept       = cell2mat({unit_basic(:).nMultiunit_kept}); 

fprintf('%s \n',subjectCode)
fprintf('nNeuron_total = %.1f, std = %.1f \n', ...
    mean(nNeuron_total), std(nNeuron_total));
fprintf('nNeuron_kept = %.1f, std = %.1f, minNum = %d, maxNum = %d \n',  ...
    mean(nNeuron_kept), std(nNeuron_kept), min(nNeuron_kept), max(nNeuron_kept));

fprintf('nSingleunit_kept = %.1f, std = %.1f, minNum = %d, maxNum = %d \n', ...
    mean(nSingleunit_kept), std(nSingleunit_kept), min(nSingleunit_kept), max(nSingleunit_kept));
fprintf('nMultiunit_kept = %.1f, std = %.1f, minNum = %d, maxNum = %d \n', ...
    mean(nMultiunit_kept), std(nMultiunit_kept), min(nMultiunit_kept), max(nMultiunit_kept));


fprintf('All units, mean SNR = %.2f, std = %.2f \n', ...
    mean(waveformSNR_neuronAll), std(waveformSNR_neuronAll))
fprintf('kept units, mean SNR = %.2f, std = %.2f \n', ...
    mean(waveformSNR_neuronKept), std(waveformSNR_neuronKept))








