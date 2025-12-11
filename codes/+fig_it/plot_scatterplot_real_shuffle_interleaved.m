function [stats_info,h] = plot_scatterplot_real_shuffle_interleaved(results_all,subjectCode, plotOptions)
global bpGlobal ftsize

tBin_toPlot = plotOptions.tBin_toPlot;
sessionType = plotOptions.sessionType; 
taskPlot    = plotOptions.taskPlot;




switch subjectCode
    case 'Ro'
        plotMarker = bpGlobal.marker.rolo;
       
    case 'Gr'
        plotMarker = bpGlobal.marker.gremlin;
end
if isfield(results_all(1),'timeWinIndex')
    idx_base =  cell2mat({results_all(:).timeWinIndex}) == tBin_toPlot & strcmp({results_all(:).sessionType},sessionType);
else
    idx_base = cell2mat({results_all(:).timeWin}) == tBin_toPlot & strcmp({results_all(:).sessionType},sessionType);

end
%%%%  cardinal only %%%%%
idx =  idx_base &  cell2mat({results_all(:).isInterleavedEpoch}) == 1 & ...
    contains({results_all(:).sessionStr},subjectCode) & strcmp({results_all(:).taskType},'cardinal');
X_cardinalOnly = cell2mat({results_all(idx).behavMeasure_cardinal});
Y_cardinalOnly_real = cell2mat({results_all(idx).combine_fisher});
Y_cardinalOnly_shuffle = cell2mat({results_all(idx).combine_shuffle});


%%%%  oblique only %%%%%
idx =  idx_base &  cell2mat({results_all(:).isInterleavedEpoch}) == 1 & ...
    contains({results_all(:).sessionStr},subjectCode) & strcmp({results_all(:).taskType},'oblique');
X_obliqueOnly = cell2mat({results_all(idx).behavMeasure_oblique});
Y_obliqueOnly_real = cell2mat({results_all(idx).combine_fisher});
Y_obliqueOnly_shuffle = cell2mat({results_all(idx).combine_shuffle});

switch taskPlot
    case 'cardinal'
        [stats_info,h] = fig.plot_two_scatter(X_cardinalOnly,Y_cardinalOnly_real,Y_cardinalOnly_shuffle,plotMarker,'Cardinal Epoch');
    case 'oblique'
        [stats_info,h] = fig.plot_two_scatter(X_obliqueOnly,Y_obliqueOnly_real,Y_obliqueOnly_shuffle,plotMarker,'Oblique Epoch');
end