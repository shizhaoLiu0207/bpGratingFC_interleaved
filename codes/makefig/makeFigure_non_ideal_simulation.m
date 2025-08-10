clear all
clc
close all

global   bpGlobal  ftsize
bpGratingFCGlobal();
save_folder = '../../figures/figures_final/model_fisher_nonideal';
%% 1. effect of delta
b_PF                = 0.00;
cardinal_delta      = 0.08;
oblique_delta       = 0.05;
cardinal_prior      = 1;
oblique_prior       = 1;
[data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);

figure
set(gcf,'Units','inches','Position',[0,0,12,4])
makefig_session_simulation(data_struct, results_fisher);
sgtitle(sprintf('b_{PF} = %d, \\color{red}Cardinal: \\delta = %.2f, prior = %.1f,  \\color{blue}Oblique: \\delta = %.2f, prior = %.1f',...
            b_PF, cardinal_delta, cardinal_prior, oblique_delta, oblique_prior),'interpreter','tex','fontweight','bold');
save_name = fullfile(save_folder,'effect_delta_simple.svg');
print(save_name,'-dsvg','-vector')
%% 2. effect of prior

b_PF                = 0.00;
cardinal_delta      = 0.08;
oblique_delta       = 0.08;
cardinal_prior      = 1;
oblique_prior       = 0.5;
[data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);

figure
set(gcf,'Units','inches','Position',[0,0,12,4])
makefig_session_simulation(data_struct, results_fisher);

sgtitle(sprintf('b_{PF} = %d, \\color{red}Cardinal: \\delta = %.2f, prior = %.1f,  \\color{blue}Oblique: \\delta = %.2f, prior = %.1f',...
            b_PF, cardinal_delta, cardinal_prior, oblique_delta, oblique_prior),'interpreter','tex','fontweight','bold');
save_name = fullfile(save_folder,'effect_prior_simple.svg');
print(save_name,'-dsvg','-vector')
%% 3. effect of b_PF
b_PF                = 0.8;
cardinal_delta      = 0.08;
oblique_delta       = 0.08;
cardinal_prior      = 1;
oblique_prior       = 1;
[data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);

figure
set(gcf,'Units','inches','Position',[0,0,12,4])
makefig_session_simulation(data_struct, results_fisher);

sgtitle(sprintf('b_{PF} = %.1f, \\color{red}Cardinal: \\delta = %.2f, prior = %.1f,  \\color{blue}Oblique: \\delta = %.2f, prior = %.1f',...
            b_PF, cardinal_delta, cardinal_prior, oblique_delta, oblique_prior),'interpreter','tex','fontweight','bold')
save_name = fullfile(save_folder,'effect_bPF_simple.svg');
print(save_name,'-dsvg','-vector')
%% 4. effect of sub_sample
b_PF                = 0;
cardinal_delta      = 0.08;
oblique_delta       = 0.08;
cardinal_prior      = 1;
oblique_prior       = 1;

session_name_str    = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior); 
session_name = ['Model_',session_name_str];

filter_name_list{1} = 'subset_32_256_random_1000';
filter_name_list{2} = 'nNeuron_256_bPF_original_0_bPF_sample_2_tao_cardinal_0_tao_oblique_0_nSample_128_random_1000.mat';
filter_name_list{3} = 'nNeuron_256_bPF_original_0_bPF_sample_-2_tao_cardinal_0_tao_oblique_0_nSample_128_random_1000.mat';
filter_name_list{4} = 'nNeuron_256_bPF_original_0_bPF_sample_0_tao_cardinal_2_tao_oblique_0_nSample_128_random_1000.mat';
%filter_name_list{5} = 'nNeuron_256_bPF_original_0_bPF_sample_0_tao_cardinal_-2_tao_oblique_0_nSample_128_random_1000.mat';
% filter_name_list{5} = 'nNeuron_256_bPF_original_0_bPF_sample_0_tao_cardinal_2_tao_oblique_2_nSample_128_random_1000.mat';
% filter_name_list{6} = 'nNeuron_256_bPF_original_0_bPF_sample_0_tao_cardinal_2_tao_oblique_-2_nSample_128_random_1000.mat';

titleStr_list{1} = 'Original';
name_str_list{1} = 'original';
for n = 2:4

    tokens = regexp(filter_name_list{n}, 'nNeuron_256_bPF_original_0_bPF_sample_([-]?[\d_]+)_tao_cardinal_([-]?[\d_]+)_tao_oblique_([-]?[\d_]+)_nSample_128_random_1000.mat', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params        = tokens{1}; % Extract matched tokens
   
    bPF_sample_str          = extracted_params{1}; 
    tao_cardinal_str        = extracted_params{2}; 
    tao_oblique_str         = extracted_params{3}; 

    titleStr_list{n} = sprintf('b_{sample} = %s  \\tau_{cardinal} = %s  \\tau_{oblique} = %s',...
                        bPF_sample_str, tao_cardinal_str, tao_oblique_str);
    name_str_list{n} = sprintf('b_%s_tao_c_%s_tao_o_%s',bPF_sample_str, tao_cardinal_str, tao_oblique_str);
end




for n = 1:4

    figure

    set(gcf,'Units','inches','Position',[0,0,12,3])
    fisher_struct = load(sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/%s/individual_sessions_cross/%s',...
            filter_name_list{n},session_name ));
    
    results_fisher = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_struct.dat_fisher_cross);
    results_fisher = get_sample_CI_cross(results_fisher);
    
    
   
    makefig_session_simulation([], results_fisher);

    sgtitle(titleStr_list{n},'fontweight','bold')

    save_name = fullfile(save_folder,sprintf('effect_sample_simple_%s.svg',name_str_list{n}));
    print(save_name,'-dsvg','-vector')
end


%% helper function
function makefig_session_simulation(data_struct, results_fisher)
 
    if ~isempty(data_struct)
        ax_1 = subplot(3,3,[1,4,7]);
        set(ax_1,'position',get(ax_1,'position')+[-0.05 0.03 0.02 -0.08]);
        plotOptions.style_cardinal = '-';
        plotOptions.style_oblique = '-';
        plotOptions.ftsize = 14;
        fig_it.plot_synth_interleaved_psycurve(data_struct.synthData_interleaved, plotOptions); 
    end
    %%%% fisher real and cross 
    if ~isempty(data_struct)
        ax_2 = subplot(3,3,[2,5]);
    else
        ax_2 = subplot(3,2,[1,3]);
    end
    set(ax_2,'position',get(ax_2,'position')+[-0.02 0.03 0.02 -0.08]);
    plotOptions = struct();
    plotOptions.errorbar = 'CI_sample';
    plotOptions.dottest = false;
    plotOptions.plotShuffle = false;
    plotOptions.ftsize = 14;
    fig_it.plot_bar_cross_Info(results_fisher, results_fisher(1).sessionStr, plotOptions); 
    
    %%%% diff cross- real
    if ~isempty(data_struct)
        ax_3 = subplot(3,3,8);
    else
        ax_3 = subplot(3,2,5);
    end
    set(ax_3,'position',get(ax_3,'position')+[-0.02 0.01 0.02 -0.03]);
    plotOptions = struct();
    plotOptions.errorbar = 'CI_sample';
    plotOptions.ftsize = 14;
    plotOptions.plotdata = 'info';
    plotOptions.markersize = 6;
    fig_it.plot_diff_errorbar(results_fisher, results_fisher(1).sessionStr, plotOptions)
    ylim([-5,60])
    
    if ~isempty(data_struct)
       ax_4 = subplot(3,3,[3,6]);
    else
       ax_4 = subplot(3,2,[2,4]);
    end
    set(ax_4,'position',get(ax_4,'position')+[0 0.03 0.02 -0.08])
    plotOptions = struct();
    plotOptions.errorbar = 'CI_sample';
    plotOptions.dottest = false;
    plotOptions.plotShuffle = false;
    plotOptions.ftsize = 14;
    plotOptions.plotPercent = true;
    fig_it.plot_bar_cross_deltaInfo(results_fisher, results_fisher(1).sessionStr, plotOptions); 
    ylim([40, 100])

    if ~isempty(data_struct)
        ax_5 = subplot(3,3,9);
    else
        ax_5 = subplot(3,2,6);
    end
    set(ax_5,'position',get(ax_5,'position')+[0 0.01 0.02 -0.03])
    plotOptions = struct();
    plotOptions.errorbar = 'CI_sample';
    plotOptions.ftsize = 14;
    plotOptions.plotdata = 'delta';
    plotOptions.markersize = 6;
    fig_it.plot_diff_errorbar(results_fisher, results_fisher(1).sessionStr, plotOptions)
    ylim([-5,25])
end

function [data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior)

    session_name_str    = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior); 
    
    data_name = ['synthData_use_interleaved_',session_name_str];
    session_name= ['Model_',session_name_str];
    
    data_folder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved/batch_1';
    fisher_folder = '../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/subset_32_256_random_1000/individual_sessions_cross';
    
    data_struct = load(fullfile(data_folder, data_name));
    fisher_struct = load(fullfile(fisher_folder, session_name));
    
    
    results_fisher = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_struct.dat_fisher_cross);
    results_fisher = get_sample_CI_cross(results_fisher);

end