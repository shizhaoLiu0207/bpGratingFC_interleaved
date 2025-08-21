function stats_info = plot_timecourse_multiple_interleaved(results_combined,subjectCode,session_list, plotOptions)
tBin_toPlot = plotOptions.tBin_toPlot;
sessionType = plotOptions.sessionType;
plotBehav   = plotOptions.plotBehav;
uncertaintyName = plotOptions.uncertaintyName;
xlimitMax = 0;
corr_type = 'spearman';
global bpGlobal ftsize
switch subjectCode
    case 'Ro'
        plotMarker = bpGlobal.marker.rolo;
        specialfill = 1;
    case 'Gr'
        plotMarker = bpGlobal.marker.gremlin;
        specialfill = 0;
end
%linearSlope(~ismember({linearSlope(:).sessionStr},session_list)) = [];
if ~isfield(results_combined(1),'timeWinIndex')
    idx_base = cell2mat({results_combined(:).timeWin}) == tBin_toPlot & strcmp({results_combined(:).sessionType},sessionType);
else
    idx_base = cell2mat({results_combined(:).timeWinIndex}) == tBin_toPlot & strcmp({results_combined(:).sessionType},sessionType);
end

%%%%%%%%%%% cardinal only %%%%%%%%%%%
idx =  idx_base &  cell2mat({results_combined(:).isInterleavedEpoch}) == 1 & ...
    contains({results_combined(:).sessionStr},subjectCode) & strcmp({results_combined(:).taskType},'cardinal');
X = cellfun(@(x)find(strcmp(session_list,x)),{results_combined(idx).sessionStr});
xlimitMax = max([xlimitMax,max(X)+1]);
%%%% cardinal real fisher
Y_cardinal          = cell2mat({results_combined(idx).combine_fisher});
Y_cardinal          = Y_cardinal';
switch uncertaintyName
    case 'std'
        Y_cardinal_error      = cell2mat({results_combined(idx).combine_fisher_std})';
    case 'CI'
        tmp = {results_combined(idx).combine_fisher_CI};
        Y_cardinal_error = cat(1,tmp{:});
end
plotOptions.taskStr = 'Cardinal';
plotOptions.epochStr = 'Cardinal';
plotOptions.plotMarker = plotMarker;
plotOptions.plotStyle = '-';
plotOptions.specialfill = 0;
plotOptions.dataStr = 'real';
plotOptions.doFill  = 1;
plotOptions.errorbarAlpha = 0.5;
plotOptions.realErrorbar = 0;
[~,h(1)] = fig.plot_time_course(X, Y_cardinal, Y_cardinal_error,plotOptions);

%%%% cardinal shuffled fisher
Y_cardinal      = cell2mat({results_combined(idx).combine_shuffle});
Y_cardinal      = Y_cardinal';
switch uncertaintyName
    case 'std'
        Y_cardinal_error      = cell2mat({results_combined(idx).combine_shuffle_std})';
    case 'CI'
        tmp = {results_combined(idx).combine_shuffle_CI};
        Y_cardinal_error = cat(1,tmp{:});
end
plotOptions.taskStr = 'Cardinal';
plotOptions.epochStr = 'Cardinal';
plotOptions.plotMarker = plotMarker;
plotOptions.plotStyle = '--';
plotOptions.specialfill = 0;
plotOptions.dataStr = 'shuffle';
plotOptions.doFill  = 0;
plotOptions.errorbarAlpha = 0.3;
plotOptions.realErrorbar = 0;
[~,h(2)] = fig.plot_time_course(X, Y_cardinal, Y_cardinal_error, plotOptions);

if plotBehav
    %%% cardinal behav fisher
    Y_cardinal      = cell2mat({results_combined(idx).behav_I})';
    switch uncertaintyName
        case 'std'
            Y_cardinal_error      = zeros(size(Y_cardinal));
        case 'CI'
            Y_cardinal_error = [zeros(size(Y_cardinal)), zeros(size(Y_cardinal))];
    end

    plotOptions.taskStr = 'Cardinal';
    plotOptions.epochStr = 'Cardinal';
    plotOptions.plotMarker = plotMarker;
    plotOptions.plotStyle = ':';
    plotOptions.specialfill = 0;
    plotOptions.dataStr = 'behav';
    plotOptions.doFill  = 0;
    plotOptions.realErrorbar = 0;
    [~,h(3)] = fig.plot_time_course(X, Y_cardinal, Y_cardinal_error, plotOptions);

  
    %  fprintf('%s, Cardinal,corr(I_real,I_behav), r = %.2f, p = %.3f \n', subjectCode, r,p);


end
[r,p] = corr(cell2mat({results_combined(idx).combine_fisher})',cell2mat({results_combined(idx).behav_I})','type',corr_type);
stats_info.r_fisher_behav_cardinal = r;
stats_info.p_fisher_behav_cardinal = p;
stats_info.df_fisher_behav_cardinal = sum(idx)-2;
%%%%%%%%%%% oblique only %%%%%%%%%%%
idx =  idx_base &  cell2mat({results_combined(:).isInterleavedEpoch}) == 1 & ...
    contains({results_combined(:).sessionStr},subjectCode) &  strcmp({results_combined(:).taskType},'oblique');
X = cellfun(@(x)find(strcmp(session_list,x)),{results_combined(idx).sessionStr});
xlimitMax = max([xlimitMax,max(X)+1]);
%%% oblique real fisher
Y_oblique = cell2mat({results_combined(idx).combine_fisher});
Y_oblique = Y_oblique';
switch uncertaintyName
    case 'std'
        Y_oblique_error      = cell2mat({results_combined(idx).combine_fisher_std})';
    case 'CI'
        tmp = {results_combined(idx).combine_fisher_CI};
        Y_oblique_error = cat(1,tmp{:});
end

plotOptions.taskStr = 'Oblique';
plotOptions.epochStr = 'Oblique';
plotOptions.plotMarker = plotMarker;
plotOptions.plotStyle = '-';
plotOptions.specialfill = specialfill;
plotOptions.dataStr = 'real';
plotOptions.doFill  = 1;
plotOptions.errorbarAlpha = 0.5;
plotOptions.realErrorbar = 0;
[~,h(4)] = fig.plot_time_course(X, Y_oblique, Y_oblique_error, plotOptions);

%%%% oblique shuffle fisher
Y_oblique = cell2mat({results_combined(idx).combine_shuffle});
Y_oblique = Y_oblique';
switch uncertaintyName
    case 'std'
        Y_oblique_error      = cell2mat({results_combined(idx).combine_shuffle_std})';
    case 'CI'
        tmp = {results_combined(idx).combine_shuffle_CI};
        Y_oblique_error = cat(1,tmp{:});
end

plotOptions.taskStr = 'Oblique';
plotOptions.epochStr = 'Oblique';
plotOptions.plotMarker = plotMarker;
plotOptions.plotStyle = '--';
plotOptions.specialfill = specialfill;
plotOptions.dataStr = 'shuffle';
plotOptions.doFill  = 0;
plotOptions.errorbarAlpha = 0.3;
plotOptions.realErrorbar = 0;
[~,h(5)] = fig.plot_time_course(X, Y_oblique, Y_oblique_error, plotOptions);

if plotBehav
    %%%%% oblique behav fisher
    Y_oblique = cell2mat({results_combined(idx).behav_I})';
    switch uncertaintyName
        case 'std'
            Y_oblique_error = zeros(size(Y_oblique));
        case 'CI'
            Y_oblique_error = [zeros(size(Y_oblique)),zeros(size(Y_oblique))];
    end

    plotOptions.taskStr = 'Oblique';
    plotOptions.epochStr = 'Oblique';
    plotOptions.plotMarker = plotMarker;
    plotOptions.plotStyle = ':';
    plotOptions.specialfill = specialfill;
    plotOptions.dataStr = 'behav';
    plotOptions.doFill  = 0;
    plotOptions.realErrorbar = 0;
    [~,h(6)] = fig.plot_time_course(X, Y_oblique, Y_oblique_error, plotOptions);

end
[r,p] = corr(cell2mat({results_combined(idx).combine_fisher})',cell2mat({results_combined(idx).behav_I})','type',corr_type);
%fprintf('%s, Oblique,corr(I_real,I_behav), r = %.2f, p = %.3f \n', subjectCode,r,p);
stats_info.r_fisher_behav_oblique = r;
stats_info.p_fisher_behav_oblique = p;
stats_info.df_fisher_behav_oblique = sum(idx)-2;

xlim([0,xlimitMax])
set(gca,'fontsize',ftsize);

xlabel('Session number');
%ylabel(ylabelStr)


box off

% if plotBehav
%     legend(h, 'Cardinal I_{real}','Cardinal I_{shuffle}','Cardinal I_{behav}','Oblique I_{real}','Oblique I_{shuffle}','Oblique I_{behav}');
% else
%     legend(h([1,2,4,5]),'Cardinal I_{real}','Cardinal I_{shuffle}','Oblique I_{real}','Oblique I_{shuffle}');
% end
ylabel('Fisher Information')


stats_info.corr_type = corr_type;

end