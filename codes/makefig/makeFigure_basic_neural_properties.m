clear all
clc
close all
%%
global   bpGlobal ftsize
ftsize = 14;
bpGratingFCGlobal();

load('../../results/neural/stimulus_dprime_twoanimals_interleaved_all_trials_coef1_hVis2_FR1.mat');
%load('../../results/neural/stimulus_dprime_twoanimals_learning_passive.mat');
%% histogram plot with the maximum coherence level
%%%% monkey R, cardinal: 1.5%, 2.5%, 4.0%, 10.0%
%%%% monkey R, oblique: 5.0%, 10.0%
%%%% monkey G, oblique: 7.5%, 15%
%%%% monkey G, cardinal: 3.5%, 7.5%, 15%
doThis  = 1;
plotSession = 'mainTask'; %% 'mainTask' or 'passiveViewing'
if doThis
    plot_field  = 'dprime_stimulus';
    save_folder = '../../figures/figures_final/psth_basic_neural';

    save_name       = fullfile(save_folder ,sprintf('%s_interleaved_histograms_highestcohr_%s.svg',plot_field, plotSession)) ;
    tex_name        = fullfile(save_folder ,sprintf('%s_interleaved_histograms_highestcohr_%s.tex',plot_field, plotSession)) ;
    
    figure
    set(gcf, 'Units','inches','Position',[0,0,11.5,3])
    switch plot_field
        case 'dprime_stimulus'
            xlabelStr = 'Stimulus $d^\prime$';
            bin_edges = [0:0.15:3];
        case 'dprime_stimulus_sign'
            xlabelStr = 'Stimulus $d^\prime$ (signed)';
            bin_edges = [-3:0.3:3];
        case 'dprime_choice'
            xlabelStr = 'Choice $d^\prime$';
            bin_edges = [0:0.1:2];
        
        case 'firingRate'
            xlabelStr = 'Firing rate';
            bin_edges = [0:4:100];

    end
    epoch_list = {'monkeyR';'monkeyG'};
    position_shift_list_top = {[-0.02 0.05 0.01 -0.04];
                           [-0.01 0.05 0.01 -0.04];
                           [0.03 0.05 0.01 -0.04];
                           [0.05 0.05 0.01 -0.04]};
    % position_shift_list_bottom = {[-0.02 0.03 0.01 -0.04];
    %                        [-0.01 0.03 0.01 -0.04];
    %                        [0.03 0.03 0.01 -0.04];
    %                        [0.05 0.03 0.01 -0.04]};
    [stats_string_ttest, stats_string_ranksum, median_mean_string, stats_string_onesample] = deal(cell(2,1));
    
    color_C = bpGlobal.color_list.color_cardinal;
    color_O = bpGlobal.color_list.color_oblique;
    color_C_light = bpGlobal.color_list.color_cardinal_light;
    color_O_light = bpGlobal.color_list.color_oblique_light;

    for k = 1:2
        epoch_name = epoch_list{k};
        switch epoch_name
            case 'monkeyR'
                session_list         = bpGlobal.rolo.session_list.switching;
                cohr_list           =   [0,3,5,10];
                epoch_string        = 'Monkey R';
            case 'monkeyG'
                session_list        = bpGlobal.gremlin.session_list.interleaved_real;
                cohr_list           = [0, 3.5,7.5,15];
                epoch_string        = 'Monkey G';
        end
       
    
    
        switch plotSession
            case 'mainTask'
                [Y_cardinal, ttest_cardinal, sessions_cardinal] = get_data(results_fprime_choice_all, session_list, plot_field, cohr_list, 'cardinal');
                [Y_oblique, ttest_oblique,sessions_oblique]   = get_data(results_fprime_choice_all, session_list, plot_field, cohr_list, 'oblique');
            case 'passiveViewing'
                % [Y_early, ttest_early] = get_data(results_fprime_choice_passive_all, session_list_early, plot_field, cohr_list_passive, plotTask);
                % [Y_late, ttest_late]   = get_data(results_fprime_choice_passive_all, session_list_late, plot_field, cohr_list_passive, plotTask);
        end
           
   
        
        
      
        
        %%%%%%  cardinal task
        ax_1 = subplot(1,4,1 + (k-1) * 2); hold on
        set(ax_1,'position',get(ax_1,'position')+position_shift_list_top{1 + (k-1) * 2});
        %data_plot  = Y_cardinal;
        h = histogram(Y_cardinal,bin_edges,'FaceColor',color_C_light,'FaceAlpha',0.2,'edgecolor',color_C_light);
        if ismember(plot_field,{'dprime_stimulus';'dprime_choice';'tuningIndex';'dprime_stimulus_sign'})
            histogram(Y_cardinal(ttest_cardinal == 1),bin_edges ,'FaceColor',color_C,'FaceAlpha',0.6,'edgecolor',color_C);
            
        end

        yl = ylim; % get current y-axis limits
        plot(median(Y_cardinal,'omitnan'), yl(2)*1.05, 'v', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor','k') % triangle for data1
        set(gca,'fontsize',14)
        %if k == 1 | k == 3
        ylabel('Num. units*','Interpreter','latex');
       % end

        xlabel(xlabelStr,'interpreter','latex')
        if ismember(plot_field,{'dprime_stimulus';'dprime_choice';'tuningIndex'})
            cardinal_stats_string = {sprintf('$Median = %.2f$',median(Y_cardinal,'omitnan'));...
                sprintf('$Fraction = %d\\%%$',round(100 * sum(ttest_cardinal) / numel(ttest_cardinal)))};
        else
            cardinal_stats_string = sprintf('$Median = %.2f$',median(Y_cardinal,'omitnan'));
        end
        axpos = get(ax_1,'position');
        w = 0.01; h = 0.01; % initial size, will auto-fit with 'FitBoxToText','on'
        x = axpos(1) + axpos(3)*0.96;
        y = axpos(2) + axpos(4)*0.98;
        hText = annotation('textbox', [x y w h], ...
            'String', cardinal_stats_string, ...
            'FitBoxToText', 'on', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontSize', ftsize,'Interpreter','latex');
    



        %%%%% oblique
        ax_2 = subplot(1,4,2 + (k-1) * 2 );  hold on
        set(ax_2,'position',get(ax_2,'position') + position_shift_list_top{2 + (k-1) * 2});
   
        histogram(Y_oblique,bin_edges,'FaceColor',color_O_light,'FaceAlpha',0.2,'edgecolor',color_O_light);
        if ismember(plot_field, {'dprime_stimulus';'dprime_choice';'tuningIndex';'dprime_stimulus_sign'})
            histogram(Y_oblique(ttest_oblique == 1),bin_edges ,'FaceColor',color_O,'FaceAlpha',0.6,'edgecolor',color_O);
        end
        yl = ylim; % get current y-axis limits
        plot(median(Y_oblique,'omitnan'), yl(2)*1.05, 'v', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k') % triangle for data1
        set(gca,'fontsize',14)
        % if k == 1 | k == 3
        %      ylabel('Num. units*','Interpreter','latex');
        % end
        xlabel(xlabelStr,'interpreter','latex')
    
         %xlabel(xlabelStr)
        if ismember(plot_field,{'dprime_stimulus';'dprime_choice';'tuningIndex'})
            oblique_stats_string = {sprintf('$Median = %.2f$',median(Y_oblique,'omitnan'));...
                sprintf('$Fraction = %d\\%%$',round(100 * sum(ttest_oblique) / numel(ttest_oblique)))};
        else
            oblique_stats_string = sprintf('$Median = %.2f$',median(Y_oblique,'omitnan'));
        end
        axpos = get(ax_2,'position');
        w = 0.01; h = 0.01; % initial size, will auto-fit with 'FitBoxToText','on'
        x = axpos(1) + axpos(3)*0.95;
        y = axpos(2) + axpos(4)*0.98;
        hText = annotation('textbox', [x y w h], ...
            'String', oblique_stats_string, ...
            'FitBoxToText', 'on', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'FontSize', ftsize ,'Interpreter','latex');
        

        %%%% statsitic test, compare between cardinal and oblique
        if ~isempty(Y_cardinal) & ~isempty(Y_oblique)
            [~,p_t, ~,stats_t] = ttest2(Y_cardinal, Y_oblique);
            [p_rank,~,stats_rank]  = ranksum(Y_cardinal, Y_oblique);
            if strcmp(plot_field, 'dprime_stimulus_sign')
                [~, p_cardinal, ~, stats_cardinal] = ttest(Y_cardinal);
                [~, p_oblique, ~, stats_oblique]   = ttest(Y_oblique);
            end
        else
            p_t = nan;
            p_rank = nan;
            stats_t.df = nan;
            stats_t.tstat = nan;
            stats_rank.ranksum = nan;
            stats_rank.zval = nan;

        end
        stats_string_ttest{k}      = sprintf('%s: $\\t(%d) = %.2f$, $p = \\num{%.2e}$; \n',...
                                            epoch_string, stats_t.df, stats_t.tstat, p_t);
        stats_string_ranksum{k}    = sprintf('%s: $W = %.2f$, $z = %.2f$, $p = \\num{%.2e}$; \n',...
                                            epoch_string, stats_rank.ranksum, stats_rank.zval, p_rank);
        if strcmp(plot_field, 'dprime_stimulus_sign')
            stats_string_onesample{k} = sprintf(['%s, cardinal one-sample t-test: $\\t(%d) = %.2f$, $p = \\num{%.2e}$; \n' ...
                                                    '%s, oblique one-sample t-test: $\\t(%d) = %.2f$, $p = \\num{%.2e}$'],...
                                                    epoch_string, stats_cardinal.df, stats_cardinal.tstat, p_cardinal,...
                                                    epoch_string, stats_oblique.df, stats_oblique.tstat, p_oblique);
        else
            stats_string_onesample{k} = '';
        end
      
        %%%%% hierarchical version of t-test that compares cardinal vs. oblique
        Y_all = [Y_cardinal;Y_oblique];
        group_all = [zeros(size(Y_cardinal));ones(size(Y_oblique))];
        session_all = [sessions_cardinal;sessions_oblique];
        lmm_results = util_it.hierarchical_ttest(Y_all, group_all, session_all);
        stats_string_lmm{k} = sprintf('%s: $\\t(%1.f) = %.2f$, $p = \\num{%.2e}$; \n',...
                                            epoch_string, lmm_results.df, lmm_results.tval, lmm_results.pval);
    end
    % 
   


    print(save_name, '-dsvg','-vector');
    %%%% generate a tex file with statistics information
    fid = fopen(tex_name,'wt');
    
    fwrite(fid, [median_mean_string{1}, median_mean_string{2}, ...
        'Independent samples t-test:',...
        stats_string_ttest{1}, stats_string_ttest{2},...
        'Wilcoxon rank sum test:',...
        stats_string_ranksum{1}, stats_string_ranksum{2}, ...
        'Onesample t-test:',...
        stats_string_onesample{1}, stats_string_onesample{2},...
        'hierarchinal lmm:',...
        stats_string_lmm{1}, stats_string_lmm{2}]);
    
    fclose(fid);
end
%% helper functions

function [Y, Y_ttest,sessions] = get_data(results_fprime_choice_all,session_list,plot_field, cohr_list, plotTask)
    [Y,Y_ttest, nOri, sessions] = deal(cell(numel(session_list),1));
    session_list_all = unique({results_fprime_choice_all(:).sessionStr});
     for n = 1:numel(session_list)
         
         switch plot_field
             case {'tuningIndex';'dprime_stimulus';'dprime_stimulus_sign'}
                idx = strcmp({results_fprime_choice_all(:).sessionStr},  session_list{n}) & ...
                    strcmp({results_fprime_choice_all(:).task}, plotTask) & ...
                    cell2mat({results_fprime_choice_all(:).cohr_level}) == max(cohr_list);
                if sum(idx>0)
                    eval(sprintf('Y{n}         = results_fprime_choice_all(idx).%s;',plot_field));
                    if strcmp(plot_field,'tuningIndex') & isfield(results_fprime_choice_all,'ttest_tuningIndex')
                        Y_ttest{n}  = results_fprime_choice_all(idx).ttest_tuningIndex;
                        nOri{n} = results_fprime_choice_all(idx).nOri;
                    else
                        Y_ttest{n}  = results_fprime_choice_all(idx).ttest_stimulus;
                    end

                end
                sessions{n} = ones(size(Y{n})) * find(strcmp(session_list_all, session_list{n}));
             case {'dprime_choice';'firingRate';'fano_per_stim'}
                 idx = strcmp({results_fprime_choice_all(:).sessionStr},  session_list{n}) & ...
                     strcmp({results_fprime_choice_all(:).task}, plotTask);
                 eval(sprintf('tmp         = {results_fprime_choice_all(idx).%s};',plot_field));
                 Y{n}     = mean(cat(2,tmp{:}), 2);
    
                 if strcmp(plot_field, 'dprime_choice')
                     tmp = {results_fprime_choice_all(idx).ttest_choice};
                     Y_ttest{n} = sum(cat(2,tmp{:}),2) > 0;
                    
                 end
                  sessions{n} = ones(size(Y{n})) * find(strcmp(session_list_all, session_list{n}));
             case {'dprime_stimulus_norm'}
             
                idx_base = strcmp({results_fprime_choice_all(:).sessionStr},  session_list{n}) ...
                            & strcmp({results_fprime_choice_all(:).task}, plotTask);
                cohr_list_real = unique(cell2mat({results_fprime_choice_all(idx_base).cohr_level}));
                cohr_list_real(cohr_list_real > 20 | cohr_list_real == 0) = []; % only use within 20% coherence level
                dprime_stimulus_norm = cell(numel(cohr_list_real),1);

                for i = 1:numel(cohr_list_real)
                    idx = idx_base & cell2mat({results_fprime_choice_all(:).cohr_level}) == cohr_list_real(i);

                    dprime_stimulus_norm{i} = results_fprime_choice_all(idx).dprime_stimulus / cohr_list_real(i);

                end
                   
                Y{n} = mean(cat(2,dprime_stimulus_norm{:}),2);
                Y_ttest{n} = nan * ones(size(Y{n}));
                sessions{n} = ones(size(Y{n})) * find(strcmp(session_list_all, session_list{n}));
         end
    
     end
    
     Y = cat(1,Y{:});
     Y_ttest = cat(1, Y_ttest{:});

     sessions = cat(1, sessions{:});
     % if exist('nOri','var')
     %    nOri = cat(1,nOri{:});
     %    Y(nOri < 8) = [];
     %    Y_ttest(nOri < 8) = [];
     % end


end
