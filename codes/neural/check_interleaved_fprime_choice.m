clear all
clc
close all
%%% This script checks stimulus tuning and choice tuning properties of
%%% neurons in the interleaved session
%%
global bpGlobal ftsize
bpGratingFCGlobal;
session_list_ro = bpGlobal.rolo.session_list.switching;

session_list_gr = bpGlobal.gremlin.session_list.interleaved_real;

%%% prepare session list and neuron filter
session_list_rolo           = bpGlobal.rolo.session_list.switching;
% session_list_rolo_good      = bpGlobal.rolo.session_list.switching_good;
% session_list_rolo_notgood   = setdiff(session_list_rolo, session_list_rolo_good);

session_list_gremlin        = bpGlobal.gremlin.session_list.interleaved_real;
% gremlin_blocksize           = bpGlobal.gremlin.session_list.interleaved_blockSize;
% session_list_gremlin_trial = gremlin_blocksize(cell2mat(gremlin_blocksize(:,2)) == 0, 1);
% session_list_gremlin_block = gremlin_blocksize(cell2mat(gremlin_blocksize(:,2)) > 0, 1);


versionName     = 'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2';
%% compute all tuning properties
session_list_all = [session_list_rolo;session_list_gremlin];
nSession = numel(session_list_all);
property_list  =  {'fprime';'dprime';'fanoFactor';'firingRate';'tuningIndex'};
results_fprime_choice_all = struct([]);
for n = 1:nSession
    fprintf('Running session %d/%d \n', n, nSession);
    sessionStr = session_list_all{n};
    data_out  = get_neural_data(sessionStr, versionName);

    results_fprime_choice = neu.get_stimulus_choice_tuning(data_out, sessionStr, property_list);
    results_fprime_choice_all = [results_fprime_choice_all,results_fprime_choice];


end
save(sprintf('../../results/neural/stimulus_dprime_twoanimals_interleaved_%s',versionName),'results_fprime_choice_all');

% %% combine across coherence
% results_fprime_choice_combined = struct([]);
% for n = 1:nSession
%     sessionStr = session_list_all{n};
%     %%%%%%%%%%%%%%%%% cardinal %%%%%%%%%%%%%%%%%%%%
%     idx = strcmp({results_fprime_choice_all(:).sessionStr}, sessionStr) & ...
%             strcmp({results_fprime_choice_all(:).task}, 'cardinal');
% 
%     %%%% fprime cardinal
%     tmp = {results_fprime_choice_all(idx).fprime};
%     fprime_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).fprime_norm};
%     fprime_norm_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%% choice signal cardinal
%     tmp = {results_fprime_choice_all(idx).choice_signal};
%     choice_signal_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).choice_signal_norm};
%     choice_signal_norm_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%% dprime cardinal
%     %%%% Shizhao on 04/27/2025: can not combine dprime this way now because  dprime was
%     %%%% unnormalized in the previous step
%     % tmp = {results_fprime_choice_all(idx).dprime_stimulus};
%     % dprime_stimulus_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).dprime_choice};
%     dprime_choice_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).dprime_stimulus_sign};
%     dprime_stimulus_sign_combine_cardinal = mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%% fano factor cardinal
%     tmp = {results_fprime_choice_all(idx).fano_per_stim};
%     fano_per_stim_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).fano_per_stim_choice};
%     fano_per_stim_choice_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%% abs(choice signal) cardinal
%     tmp = {results_fprime_choice_all(idx).abs_choice_signal};
%     abs_choice_signal_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).abs_choice_signal_norm};
%     abs_choice_signal_norm_combine_cardinal =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%%%%%%%%%%%%%%% oblique %%%%%%%%%%%%%%%%%%%%
%     idx = strcmp({results_fprime_choice_all(:).sessionStr}, sessionStr) & ...
%             strcmp({results_fprime_choice_all(:).task}, 'oblique');
% 
%     %%%% fprime oblique
%     tmp = {results_fprime_choice_all(idx).fprime};
%     fprime_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).fprime_norm};
%     fprime_norm_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%% choice signal oblique
%     tmp = {results_fprime_choice_all(idx).choice_signal};
%     choice_signal_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).choice_signal_norm};
%     choice_signal_norm_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%% dprime oblique
%     %%%% Shizhao on 04/27/2025: can not combine dprime this way now because  dprime was
%     %%%% unnormalized in the previous step
%     % tmp = {results_fprime_choice_all(idx).dprime_stimulus};
%     % dprime_stimulus_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).dprime_choice};
%     dprime_choice_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).dprime_stimulus_sign};
%     dprime_stimulus_sign_combine_oblique = mean(cat(2, tmp{:}),2 ,'omitnan');
%     %%%% fano factor oblique
%     tmp = {results_fprime_choice_all(idx).fano_per_stim};
%     fano_per_stim_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).fano_per_stim_choice};
%     fano_per_stim_choice_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
%     %%%% abs(choice signal) oblique
%     tmp = {results_fprime_choice_all(idx).abs_choice_signal};
%     abs_choice_signal_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
%     tmp = {results_fprime_choice_all(idx).abs_choice_signal_norm};
%     abs_choice_signal_norm_combine_oblique =  mean(cat(2, tmp{:}),2 ,'omitnan');
% 
% 
%     results_fprime_choice_combined(n).sessionStr                            = sessionStr;
%     results_fprime_choice_combined(n).fprime_combine_cardinal               = fprime_combine_cardinal;
%     results_fprime_choice_combined(n).fprime_combine_oblique                = fprime_combine_oblique;
%     results_fprime_choice_combined(n).fprime_norm_combine_cardinal          = fprime_norm_combine_cardinal;
%     results_fprime_choice_combined(n).fprime_norm_combine_oblique           = fprime_norm_combine_oblique;
%     results_fprime_choice_combined(n).choice_signal_combine_cardinal        = choice_signal_combine_cardinal;
%     results_fprime_choice_combined(n).choice_signal_combine_oblique         = choice_signal_combine_oblique;
%     results_fprime_choice_combined(n).choice_signal_norm_combine_cardinal   = choice_signal_norm_combine_cardinal;
%     results_fprime_choice_combined(n).choice_signal_norm_combine_oblique    = choice_signal_norm_combine_oblique;
% 
%     results_fprime_choice_combined(n).abs_choice_signal_combine_cardinal        = abs_choice_signal_combine_cardinal;
%     results_fprime_choice_combined(n).abs_choice_signal_combine_oblique         = abs_choice_signal_combine_oblique;
%     results_fprime_choice_combined(n).abs_choice_signal_norm_combine_cardinal   = abs_choice_signal_norm_combine_cardinal;
%     results_fprime_choice_combined(n).abs_choice_signal_norm_combine_oblique    = abs_choice_signal_norm_combine_oblique;
% 
%     results_fprime_choice_combined(n).dprime_stimulus_combine_cardinal      = dprime_stimulus_combine_cardinal;
%     results_fprime_choice_combined(n).dprime_stimulus_combine_oblique       = dprime_stimulus_combine_oblique;
%     results_fprime_choice_combined(n).dprime_choice_combine_cardinal        = dprime_choice_combine_cardinal;
%     results_fprime_choice_combined(n).dprime_choice_combine_oblique         = dprime_choice_combine_oblique;
%     results_fprime_choice_combined(n).dprime_stimulus_sign_combine_cardinal = dprime_stimulus_sign_combine_cardinal;
%     results_fprime_choice_combined(n).dprime_stimulus_sign_combine_oblique  = dprime_stimulus_sign_combine_oblique;
%     results_fprime_choice_combined(n).fano_per_stim_combine_cardinal        = fano_per_stim_combine_cardinal;
%     results_fprime_choice_combined(n).fano_per_stim_combine_oblique         = fano_per_stim_combine_oblique;
%     results_fprime_choice_combined(n).fano_per_stim_choice_combine_cardinal = fano_per_stim_choice_combine_cardinal;
%     results_fprime_choice_combined(n).fano_per_stim_choice_combine_oblique  = fano_per_stim_choice_combine_oblique;
% 
% end
% %% make figures
% doPlot = 0;
% saveFolder = '../../figures/figures_informal/choice_fprime_fano_interleaved';
% subject_list  = {'MonkeyG'};
% plotData_list = {'fprime'; 'fprime_norm'; 'choice_signal'; 'choice_signal_norm';
%                'abs_choice_signal';'abs_choice_signal_norm'
%                 'dprime_stimulus' ; 'dprime_choice';
%                'fano_per_stim_choice';
%                  };
% 
% for i = 1:numel(subject_list)
%     for j = 1:numel(plotData_list)
%         subject =  subject_list{i};
%         plotData = plotData_list{j};
%         switch subject
%             case 'MonkeyR'
%                 session_list_plot = session_list_rolo;
%             case 'MonkeyG'
%                 session_list_plot = session_list_gremlin_trial;
%         end
% 
%         idx = ismember({results_fprime_choice_combined(:).sessionStr}, session_list_plot);
%         eval(sprintf('tmp = {results_fprime_choice_combined(idx).%s_combine_cardinal};',plotData));
% 
%         data_cardinal = cat(1, tmp{:});
% 
%         eval(sprintf('tmp = {results_fprime_choice_combined(idx).%s_combine_oblique};',plotData));
%         data_oblique = cat(1,tmp{:});
% 
% 
% 
%         t_test_str = plot_scatter_histogram(data_cardinal, data_oblique, bpGlobal.color_list.color_cardinal, bpGlobal.color_list.color_oblique, doPlot);
%         if doPlot
%             figure;
%             set(gcf,'unit','inch','position',[0,0,7,5])
%             ftsize = 16;
%             plot_scatter_histogram(data_cardinal, data_oblique, bpGlobal.color_list.color_cardinal, bpGlobal.color_list.color_oblique, doPlot);
%             sgtitle(sprintf('%s, %s', subject, plotData),'interpreter','none','fontsize', ftsize + 2,'fontweight','bold');
%             saveName = fullfile(saveFolder, sprintf('%s_%s',subject,plotData));
%             print(saveName,'-dsvg','-vector')
%         else
%             t_test_str = plot_scatter_histogram(data_cardinal, data_oblique, bpGlobal.color_list.color_cardinal, bpGlobal.color_list.color_oblique, doPlot);
%             fprintf('%s_%s: %s.\n\n',subject,plotData,t_test_str);
% 
%         end
%     end
% end
% %% more detailed figures on dprime_stimulus
% 
% plotOption = 'abs'; % abs or signed
% 
% figure;
% %%%% plot monkey R
% idx = ismember({results_fprime_choice_combined(:).sessionStr}, session_list_rolo);
% switch plotOption
%     case 'abs'
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_cardinal};
%         dprime_stimulus_cardinal = cat(1,tmp{:});
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_oblique};
%         dprime_stimulus_oblique = cat(1,tmp{:});
%     case 'signed'
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_cardinal};
%         dprime_stimulus_cardinal = cat(1,tmp{:});
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_oblique};
%         dprime_stimulus_oblique = cat(1,tmp{:});
% end
% 
% subplot(2,2,1)
% h_c = histogram(dprime_stimulus_cardinal,'facecolor','red','Normalization','probability');
% hold on
% histogram(dprime_stimulus_oblique,h_c.BinEdges,'facecolor','blue','Normalization','probability');
% box off
% xlabel('dprime stimulus');ylabel('Norm. Probability')
% [~,p,~,stats] = ttest2(dprime_stimulus_cardinal, dprime_stimulus_oblique);
% index_stimbias_mean = (mean(dprime_stimulus_cardinal) - mean(dprime_stimulus_oblique)) / ...
%                        (mean(dprime_stimulus_cardinal) + mean(dprime_stimulus_oblique)); 
% title({'Monkey R';sprintf('t-test = %.2f, p = %.1e',stats.tstat,p);...
%         sprintf('Index stimBias (mean) = %.2f',index_stimbias_mean)});
% set(gca,'fontsize',16)
% 
% subplot(2,2,2)
% h_c = cdfplot(dprime_stimulus_cardinal); hold on
% h_o = cdfplot(dprime_stimulus_oblique); 
% h_c.Color = 'red'; h_c.LineWidth = 2;
% h_o.Color = 'blue'; h_o.LineWidth = 2;
% [p,~,stats] = ranksum(dprime_stimulus_cardinal, dprime_stimulus_oblique);
% index_stimbias_median = (median(dprime_stimulus_cardinal) - median(dprime_stimulus_oblique)) / ...
%                        (median(dprime_stimulus_cardinal) + median(dprime_stimulus_oblique)); 
% title({'Monkey R';sprintf('t-test = %.2f, p = %.1e',stats.zval,p);...
%         sprintf('Index stimBias (median) = %.2f',index_stimbias_median)});
% box off
% xlabel('dprime stimulus')
% set(gca,'fontsize',16)
% %%%% plot monkey G
% idx = ismember({results_fprime_choice_combined(:).sessionStr}, session_list_gremlin);
% switch plotOption
%     case 'abs'
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_cardinal};
%         dprime_stimulus_cardinal = cat(1,tmp{:});
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_oblique};
%         dprime_stimulus_oblique = cat(1,tmp{:});
%     case 'signed'
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_cardinal};
%         dprime_stimulus_cardinal = cat(1,tmp{:});
%         tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_oblique};
%         dprime_stimulus_oblique = cat(1,tmp{:});
% end
% 
% subplot(2,2,3)
% h_c = histogram(dprime_stimulus_cardinal, 'facecolor','red','Normalization','probability');;
% hold on
% histogram(dprime_stimulus_oblique, h_c.BinEdges, 'facecolor','blue','Normalization','probability');
% box off
% xlabel('dprime stimulus');ylabel('Norm. Probability')
% [~,p,~,stats] = ttest2(dprime_stimulus_cardinal, dprime_stimulus_oblique);
% index_stimbias_mean = (mean(dprime_stimulus_cardinal) - mean(dprime_stimulus_oblique)) / ...
%                        (mean(dprime_stimulus_cardinal) + mean(dprime_stimulus_oblique)); 
% title({'Monkey G';sprintf('t-test = %.2f, p = %.1e',stats.tstat,p);...
%         sprintf('Index stimBias (mean) = %.2f',index_stimbias_mean)});
% set(gca,'fontsize',16)
% 
% subplot(2,2,4)
% h_c = cdfplot(dprime_stimulus_cardinal); hold on
% h_o = cdfplot(dprime_stimulus_oblique); 
% h_c.Color = 'red'; h_c.LineWidth = 2;
% h_o.Color = 'blue'; h_o.LineWidth = 2;
% [p,~,stats] = ranksum(dprime_stimulus_cardinal, dprime_stimulus_oblique);
% index_stimbias_median = (median(dprime_stimulus_cardinal) - median(dprime_stimulus_oblique)) / ...
%                        (median(dprime_stimulus_cardinal) + median(dprime_stimulus_oblique)); 
% title({'Monkey G';sprintf('t-test = %.2f, p = %.1e',stats.zval,p);...
%         sprintf('Index stimBias (median) = %.2f',index_stimbias_median)});
% box off
% xlabel('dprime stimulus')
% set(gca,'fontsize',16)
% 
% switch plotOption
%     case 'abs'
%        sgtitleStr = 'absolute dprime';
%     case 'signed'
%        sgtitleStr = 'signed dprime';
% end
% sgtitle(sgtitleStr,'fontsize',20);
% %% save dprime variable for future use
% idx = ismember({results_fprime_choice_combined(:).sessionStr}, session_list_rolo);
% 
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_cardinal};
% dprime_stimulus_cardinal_rolo = cat(1,tmp{:});
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_oblique};
% dprime_stimulus_oblique_rolo = cat(1,tmp{:});
% 
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_cardinal};
% dprime_stimulus_sign_cardinal_rolo = cat(1,tmp{:});
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_oblique};
% dprime_stimulus_sign_oblique_rolo = cat(1,tmp{:});
% 
% idx = ismember({results_fprime_choice_combined(:).sessionStr}, session_list_gremlin);
% 
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_cardinal};
% dprime_stimulus_cardinal_gremlin = cat(1,tmp{:});
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_combine_oblique};
% dprime_stimulus_oblique_gremlin = cat(1,tmp{:});
% 
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_cardinal};
% dprime_stimulus_sign_cardinal_gremlin = cat(1,tmp{:});
% tmp = {results_fprime_choice_combined(idx).dprime_stimulus_sign_combine_oblique};
% dprime_stimulus_sign_oblique_gremlin = cat(1,tmp{:});
% 
% saveName = '../../results/neural/stimulus_dprime_twoanimals_interleaved';
% save(saveName,'dprime_stimulus_cardinal_rolo','dprime_stimulus_oblique_rolo',...
%                'dprime_stimulus_sign_cardinal_rolo','dprime_stimulus_sign_oblique_rolo',...
%                'dprime_stimulus_cardinal_gremlin','dprime_stimulus_oblique_gremlin',...
%                'dprime_stimulus_sign_cardinal_gremlin','dprime_stimulus_sign_oblique_gremlin')
%%
function title_text = plot_scatter_histogram(x, y, x_color, y_color, doPlot)
    global ftsize
    idx_nan = isnan(x) | isnan(y);
    x(idx_nan) = [];
    y(idx_nan) = [];
    % Perform t-test
    [~, p, ~, stats] = ttest(x, y);
   
    title_text = sprintf('t-test: t = %.3f, p = %.3f', stats.tstat, p);
    if doPlot
        % Create figure with subplots
        
        
        % Define subplot layout
        scatter_ax = subplot(3,3,[2,3,5,6]); % Scatter plot
        hist_x_ax = subplot(3,3,[8,9]); % Histogram of x
        hist_y_ax = subplot(3,3,[1,4]); % Histogram of y
        
        % Compute common axis limits
        x_min = min([x; y]);
        x_max = max([x; y]);
        x_lim = [x_min, x_max];
        
        
        % Scatter plot
        axes(scatter_ax);
        scatter(x, y, 'filled'); hold on;
        plot([min([x;y]), max([x;y])], [min([x;y]), max([x;y])], 'k--'); % Diagonal line
        % xlabel('');
        % ylabel('Oblique');
        lgd =  legend(title_text);
        set(lgd,'Box','off');
        xlim(x_lim);
        ylim(x_lim); % Match y limits with x limits
        set(gca, 'fontsize', ftsize)
        box off
    
        % Histogram of x on bottom
        axes(hist_x_ax);
        histogram(x, 'FaceColor', x_color);
        lgd = legend(sprintf('Mean = %.2f', mean(x,'omitnan')));
        set(lgd,'Box','off');
        xlabel('Cardinal');
        ylabel('Count');
        xlim(x_lim); % Match x-lim with scatter plot
        set(gca, 'fontsize', ftsize)
        box off
        % Histogram of y on left
        axes(hist_y_ax);
        histogram(y, 'FaceColor', y_color);
        lgd = legend(sprintf('Mean = %.2f', mean(y,'omitnan')));
        set(lgd,'Box','off');
        xlabel('Oblique');
        ylabel('Count');
        xlim(x_lim); % Match y-lim with scatter plot
        view([90 -90]); % Rotate histogram to match y-axis
        set(gca, 'fontsize', ftsize)
        box off
    
        % % Adjust layout
        % set(scatter_ax, 'Position', [0.3, 0.3, 0.6, 0.6]);
        % set(hist_x_ax, 'Position', [0.3, 0.1, 0.6, 0.2]);
        % set(hist_y_ax, 'Position', [0.1, 0.3, 0.2, 0.6]);
        
        hold off;
    end

end
