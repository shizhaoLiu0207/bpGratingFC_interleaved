clear all
clc
close all
%%
global bpGlobal
bpGratingFCGlobal;
session_list_ro = bpGlobal.rolo.session_list.switching;

session_list_gr = bpGlobal.gremlin.session_list.interleaved_real;

%%% prepare session list and neuron filter
session_list_rolo           = bpGlobal.rolo.session_list.switching;
session_list_rolo_good      = bpGlobal.rolo.session_list.switching_good;
session_list_rolo_notgood   = setdiff(session_list_rolo, session_list_rolo_good);

session_list_gremlin        = bpGlobal.gremlin.session_list.interleaved_real;
gremlin_blocksize           = bpGlobal.gremlin.session_list.interleaved_blockSize;
session_list_gremlin_trial = gremlin_blocksize(cell2mat(gremlin_blocksize(:,2)) == 0, 1);
session_list_gremlin_block = gremlin_blocksize(cell2mat(gremlin_blocksize(:,2)) > 0, 1);


filter_name     = 'all_trials_coef1_hVis2_FR1';
%% plot psth of many different conditions
%%%
save_folder = '../../figures/figures_final/psth_basic_neural';
plot_epoch_list = {'monkeyR_all';'monkeyG_all';'monkeyR_good';'monkeyG_block';'monkeyR_not_good';'monkeyG_trial'};
for  n = 1:numel(plot_epoch_list)
    plot_epoch = plot_epoch_list{n};
    save_name = fullfile(save_folder,sprintf('psth_interleaved_%s',plot_epoch));
    switch plot_epoch
        case 'monkeyR_all'
            session_list_plot = session_list_rolo;
            sgtitleStr = 'Monkey R, all sessions';
        case 'monkeyG_all'
            session_list_plot = session_list_gremlin;
            sgtitleStr = 'Monkey G, all sessions';
        case 'monkeyR_good'
            session_list_plot = session_list_rolo_good;
            sgtitleStr = 'Monkey R, good sessions';
        case 'monkeyG_block'
            session_list_plot = session_list_gremlin_block;
            sgtitleStr = 'Monkey G, blockwise sessions';
        case 'monkeyR_not_good'
            session_list_plot = session_list_rolo_notgood;
            sgtitleStr = 'Monkey R, not-good sessions';
        case 'monkeyG_trial'
            session_list_plot = session_list_gremlin_trial;
            sgtitleStr = 'Monkey G, trial-by-trial sessions';
    end
    figure
    ftsize = 14;
    set(gcf,'unit','inch','position',[0,0,12,3.5])
    plotOptions      = struct();

    plotOptions.alpha = 0.05;
    plotOptions.n_permutations = 1000;
    plotOptions.doAutoPlot = 0;

    ax_1 = subplot(1,4,1);
    set(ax_1,'position',get(ax_1,'position')+[-0.04 0.05 0.01 -0.1]);
    plotOptions.task = 'cardinal';
    plotOptions.cohr = 'nonzero';
    fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    title('Cardinal nonzero','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');

    ax_2 = subplot(1,4,2);
    set(ax_2,'position',get(ax_2,'position')+[-0.03 0.05 0.01 -0.1]);
    plotOptions.task = 'oblique';
    plotOptions.cohr = 'nonzero';
    fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    title('Oblique nonzero','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');


    ax_3 = subplot(1,4,3);
    set(ax_3,'position',get(ax_3,'position')+[0.01 0.05 0.01 -0.1]);
    plotOptions.task = 'cardinal';
    plotOptions.cohr = 'zero';
    fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    title('Cardinal zero','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');


    ax_4 = subplot(1,4,4);
    set(ax_4,'position',get(ax_4,'position')+[0.03 0.05 0.01 -0.1]);
    plotOptions.task = 'oblique';
    plotOptions.cohr = 'zero';
    fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    title('Oblique zero','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');


    sgtitle(sgtitleStr,'fontsize',16,'fontWeight','bold')
    print(save_name,'-dsvg','-vector')
end
%% plot less conditions for paper figure
save_folder = '../../figures/figures_final/psth_basic_neural';
plot_epoch_list = {'monkeyR_all';'monkeyG_all'};
plotNorm = true;
for  n = 1:numel(plot_epoch_list)
    plot_epoch = plot_epoch_list{n};
    if plotNorm
        save_name = fullfile(save_folder,sprintf('psth_interleaved_%s_normalized',plot_epoch));
    else
        save_name = fullfile(save_folder,sprintf('psth_interleaved_%s',plot_epoch));
    end
    switch plot_epoch
        case 'monkeyR_all'
            session_list_plot = session_list_rolo;
            sgtitleStr = 'Monkey R';
        case 'monkeyG_all'
            session_list_plot = session_list_gremlin;
            sgtitleStr = 'Monkey G';
    end
    figure
    ftsize = 14;
    set(gcf,'unit','inch','position',[0,0,6,4])
    plotOptions      = struct();

    plotOptions.t    = [-50:1700]; 
    plotOptions.cohr = 'nonzero'; 
    plotOptions.plotNorm = plotNorm;
    
    plotOptions.doStats         = true; 
    plotOptions.alpha           = 0.05;
    plotOptions.n_permutations  = 1000;
    plotOptions.doAutoPlot      = 0;

    ax_1 = subplot(1,2,1);
    set(ax_1,'position',get(ax_1,'position')+[-0.04 0.05 0.01 -0.1]);
    plotOptions.task = 'cardinal';

    fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    title('Cardinal nonzero','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');

    ax_2 = subplot(1,2,2);
    set(ax_2,'position',get(ax_2,'position')+[-0.03 0.05 0.01 -0.1]);
    plotOptions.task = 'oblique';

    fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    title('Oblique nonzero','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');


    % ax_3 = subplot(1,4,3);
    % set(ax_3,'position',get(ax_3,'position')+[0.01 0.05 0.01 -0.1]);
    % plotOptions.task = 'cardinal';
    % plotOptions.cohr = 'zero';
    % fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    % title('Cardinal zero','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
    % 
    % 
    % ax_4 = subplot(1,4,4);
    % set(ax_4,'position',get(ax_4,'position')+[0.03 0.05 0.01 -0.1]);
    % plotOptions.task = 'oblique';
    % plotOptions.cohr = 'zero';
    % fig.plot_psth_norm_grouped_new(filter_name, session_list_plot, plotOptions);
    % title('Oblique zero','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');


    sgtitle(sgtitleStr,'fontsize',16,'fontWeight','bold')
    print(save_name,'-dsvg','-vector')
end