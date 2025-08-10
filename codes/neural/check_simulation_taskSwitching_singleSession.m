clear all
clc
%close all

saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved';
fileName_list{1} = 'synthData_use_interleaved_bPF_0_00_cardinal_delta_0_08_prior_1_00_oblique_delta_0_04_prior_1_00';
fileName_list{2} = 'synthData_use_interleaved_bPF_0_00_cardinal_delta_0_08_prior_1_00_oblique_delta_0_08_prior_0_50';
fileName_list{3} = 'synthData_use_interleaved_bPF_0_80_cardinal_delta_0_08_prior_1_00_oblique_delta_0_08_prior_1_00';



plotOptions.style_cardinal = '-';
plotOptions.style_oblique = '-';


figure      
set(gcf,'unit','inches','position',[0,0,15,4])
save_folder = '../../figures/figures_informal/cross_fisher';
for i = 1:3
    load(fullfile(saveFolder, fileName_list{i}));
    subplot(1,3,i)
    plot_interleaved_psycurve(synthData_interleaved, plotOptions);
    set(gca,'fontsize',20);
    box off
    % if i == 1
    %     ylabel('Percent of 0^\degree/45^\degree choice')
    % end
    xlabel('Signal level (%)')
end
save_name = fullfile(save_folder,'example_newsimulation_psycurve.svg');
print(save_name, '-dsvg');
function h = plot_interleaved_psycurve(synthData_interleaved, plotOptions)

idx_cardinal                = strcmp({synthData_interleaved(:).image_task}, 'cardinal');
decision_cardinal           = {synthData_interleaved(idx_cardinal).decision};
decision_cardinal           = cat(2,decision_cardinal{:}); 
contrast_signed_cardinal    = cell2mat({synthData_interleaved(idx_cardinal).contrast_signed});

idx_oblique                 = strcmp({synthData_interleaved(:).image_task}, 'oblique');
decision_oblique            = {synthData_interleaved(idx_oblique).decision};
decision_oblique            = cat(2,decision_oblique{:}); 
contrast_signed_oblique     = cell2mat({synthData_interleaved(idx_oblique).contrast_signed});

probChoice2_cardinal        = sum(decision_cardinal == 2, 1) / size(decision_cardinal,1);
semChoice_ori2_cardinal     = sqrt(probChoice2_cardinal .* (1 - probChoice2_cardinal) / size(decision_cardinal,1));
probChoice2_oblique         = sum(decision_oblique == 2, 1)   / size(decision_oblique,1);
semChoice_ori2_oblique      = sqrt(probChoice2_oblique .* (1 - probChoice2_oblique) / size(decision_oblique,1));


[~, i_c] = sort(contrast_signed_cardinal ,'ascend');
[~, i_o] = sort(contrast_signed_oblique ,'ascend');

h(1) = errorbar(contrast_signed_cardinal(i_c), probChoice2_cardinal(i_c), semChoice_ori2_cardinal(i_c),'LineWidth',2,'color','red','LineStyle',plotOptions.style_cardinal); hold on
h(2) = errorbar(contrast_signed_oblique(i_o),  probChoice2_oblique(i_o), semChoice_ori2_oblique(i_o),'LineWidth',2,'color','blue','LineStyle',plotOptions.style_oblique); hold on


end