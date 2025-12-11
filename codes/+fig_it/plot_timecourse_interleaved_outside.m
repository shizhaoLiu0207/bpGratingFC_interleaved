function plot_timecourse_interleaved_outside(results_combined,subjectCode,session_list,plotOptions)


tBin_toPlot = plotOptions.tBin_toPlot;
sessionType = plotOptions.sessionType;
dataName    = plotOptions.dataName;
uncertaintyName = plotOptions.uncertaintyName;
plotInterleaved = plotOptions.plotInterleaved;

xlimitMax = 0;
global ftsize bpGlobal
%tBin_toPlot = 0;
switch subjectCode
    case 'Ro'
        plotMarker = bpGlobal.marker.rolo;
        specialfill = 1;
    case 'Gr'
        plotMarker = bpGlobal.marker.gremlin;
        specialfill = 0;
end
%linearSlope(~ismember({linearSlope(:).sessionStr},session_list)) = [];
if isfield(results_combined(1),'timeWinIndex')
    idx_base = cell2mat({results_combined(:).timeWinIndex}) == tBin_toPlot & strcmp({results_combined(:).sessionType},sessionType);
else
    idx_base = cell2mat({results_combined(:).timeWin}) == tBin_toPlot & strcmp({results_combined(:).sessionType},sessionType);
end

%%%% cardinal only %%%%%
idx =  idx_base &  cell2mat({results_combined(:).isInterleavedEpoch}) == 1 & ...
    contains({results_combined(:).sessionStr},subjectCode) & strcmp({results_combined(:).taskType},'cardinal');

eval(sprintf('Y_cardinal = cell2mat({results_combined(idx).%s});', dataName));
Y_cardinal = Y_cardinal';
eval(sprintf('tmp = {results_combined(idx).%s_%s};',dataName, uncertaintyName));
Y_cardinal_errorbar = cat(1,tmp{:});
%eval(sprintf('Y_cardinal_errorbar = cell2mat({results_combined(idx).%s_%s});', dataName, uncertaintyName))

X = cellfun(@(x)find(strcmp(session_list,x)),{results_combined(idx).sessionStr});
xlimitMax = max([xlimitMax,max(X)+1]);

plotOptions.taskStr = 'Cardinal';
plotOptions.epochStr = 'Cardinal';
plotOptions.plotMarker = plotMarker;
plotOptions.plotStyle = '-';
plotOptions.specialfill = 0;
plotOptions.dataStr = 'delta';
plotOptions.doFill  = 1;
plotOptions.realErrorbar = 1;
plotOptions.errorbarAlpha = 0.4;
h(1) = fig.plot_time_course(X, Y_cardinal, Y_cardinal_errorbar, plotOptions);
%%% oblique only main task
idx =  idx_base &  cell2mat({results_combined(:).isInterleavedEpoch}) == 1 & ...
    contains({results_combined(:).sessionStr},subjectCode) & strcmp({results_combined(:).taskType},'oblique');


eval(sprintf('Y_oblique = cell2mat({results_combined(idx).%s});', dataName));
Y_oblique = Y_oblique';
eval(sprintf('tmp = {results_combined(idx).%s_%s};',dataName, uncertaintyName));
Y_oblique_errorbar = cat(1,tmp{:});
%eval(sprintf('Y_oblique_errorbar = cell2mat({results_combined(idx).%s_std});', dataName))

X = cellfun(@(x)find(strcmp(session_list,x)),{results_combined(idx).sessionStr});

plotOptions.taskStr = 'Oblique';
plotOptions.epochStr = 'Oblique';
plotOptions.plotMarker = plotMarker;
plotOptions.plotStyle = '-';
plotOptions.specialfill = specialfill;
plotOptions.dataStr = 'delta';
plotOptions.doFill  = 1;
plotOptions.realErrorbar = 1;
h(2) = fig.plot_time_course(X, Y_oblique, Y_oblique_errorbar,plotOptions);
xlimitMax = max([xlimitMax,max(X)+1]);




set(gca,'fontsize',ftsize);

xlabel('Session number');
%ylabel(ylabelStr)
xlim([0,xlimitMax])
%ylim([-0.5,1.5])
box off

%legend(h, 'Cardinal task','Oblique task');

% slope_delta_norm
% slope_fisher_norm
% slope_shuffle_norm
% slope_unit_norm
switch dataName
    case 'slope_delta_norm'
        ylabel('\Delta Fisher Information (Norm.)')
    case 'slope_fisher_norm'
         ylabel('Real Fisher Information (Norm.)')
    case 'slope_shuffle_norm'
         ylabel('Shuffled Fisher Information (Norm.)')
    case 'slope_unit_norm'
         ylabel('Individual unit Fisher Information (Norm.)')
end
end