clear all
clc
close all
%%%% This scripts plots psychometric curves of each interleaved session for
%%%% large-scale simulation
%%
saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved';
filename_list = dir(fullfile(saveFolder, '*.mat'));
nFile = numel(filename_list);
[bPF, cardinal_delta, cardinal_prior, oblique_delta, oblique_prior] = deal(zeros(nFile,1));
[bPF_str, cardinal_delta_str, cardinal_prior_str, oblique_delta_str, oblique_prior_str] = deal(cell(nFile, 1));
for n = 1:nFile
% Extract numbers using regular expressions
    tokens = regexp(filename_list(n).name, 'synthData_use_interleaved_bPF_([-]?[\d_]+)_cardinal_delta_([\d_]+)_prior_([\d_]+)_oblique_delta_([\d_]+)_prior_([\d_]+)', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params        = tokens{1}; % Extract matched tokens
    
    bPF_str{n}                 = strrep(extracted_params{1}, '_', '.'); % Replace _ with .
    cardinal_delta_str{n}      = strrep(extracted_params{2}, '_', '.');
    cardinal_prior_str{n}      = strrep(extracted_params{3}, '_', '.');
    oblique_delta_str{n}       = strrep(extracted_params{4}, '_', '.');
    oblique_prior_str{n}       = strrep(extracted_params{5}, '_', '.');


    
    bPF(n)                  = str2double(bPF_str{n});
    cardinal_delta(n)       = str2double(cardinal_delta_str{n});
    cardinal_prior(n)       = str2double(cardinal_prior_str{n});
    oblique_delta(n)        = str2double(oblique_delta_str{n});
    oblique_prior(n)        = str2double(oblique_prior_str{n});
end
[~,i_sort] = sort(bPF,'ascend');
bPF = bPF(i_sort);
cardinal_delta = cardinal_delta(i_sort);
cardinal_prior = cardinal_prior(i_sort);
oblique_delta = oblique_delta(i_sort);
oblique_prior = oblique_prior(i_sort);
filename_list = filename_list(i_sort);
%% psychometric curves of each session with different parameter settings
figFolder = '../../figures/figures_informal/Bayesian_model_simulation';
bPF_list = unique(bPF); 
prior_list = unique(cardinal_prior);
delta_list = unique(cardinal_delta);
nPrior = numel(prior_list);
for  i = 1:numel(bPF_list)
    figure
    set(gcf,'units','normalized','position',[0,0,1,1]);
    idx_base = bPF == bPF_list(i);
    for j1 = 1:nPrior
        for j2 = 1:nPrior
            subplot(nPrior, nPrior, (j1-1) * nPrior + j2)
            idx  = idx_base & cardinal_prior == prior_list(j1) & oblique_prior == prior_list(j2);
            idx_list  = find(idx);
            for n = 1:numel(idx_list)
                load(fullfile(saveFolder, filename_list(idx_list(n)).name));
                delta_c_plot = cardinal_delta(idx_list(n));
                delta_o_plot = oblique_delta(idx_list(n));
                if delta_c_plot == delta_list(1)
                    plotOptions.style_cardinal = '--';
                else
                    plotOptions.style_cardinal = '-';
                end
                if delta_o_plot == delta_list(1)
                    plotOptions.style_oblique = '--';
                else
                    plotOptions.style_oblique = '-';
                end
                 plot_interleaved_psycurve(synthData_interleaved, plotOptions);
                 
            end
            set(gca,'fontsize',16)
            box off
            title(sprintf('Cardinal prior = %.2f, Oblique prior = %.2f',prior_list(j1), prior_list(j2)))
           
        end
    end
    sgtitle(sprintf('b PF = %.1f',bPF_list(i)),'fontsize',20);
    savename = fullfile(figFolder, sprintf('psyCurve_real_nonuniform_bPF_%s.png',bPF_str{idx_list(n)}));
    print(gcf, savename,'-dpng')
    

end
%% Check accuracy at contrast = 6 with different parameters
[accuracy_cardinal, accuracy_oblique, semAcc_cardinal, semAcc_oblique] = deal(zeros(numel(filename_list),1));
for n = 1:numel(filename_list)
    load(fullfile(saveFolder, filename_list(n).name));

    nCorrect_cardinal = 0; nTotal_cardinal = 0;
    idx = strcmp({synthData_interleaved(:).image_task},'cardinal') & cell2mat({synthData_interleaved(:).contrast_signed}) == 6;
    nCorrect_cardinal = nCorrect_cardinal + sum(synthData_interleaved(idx).decision == 2);
    nTotal_cardinal   = nTotal_cardinal + numel(synthData_interleaved(idx).decision);
    idx = strcmp({synthData_interleaved(:).image_task},'cardinal') & cell2mat({synthData_interleaved(:).contrast_signed}) == -6;
    nCorrect_cardinal = nCorrect_cardinal + sum(synthData_interleaved(idx).decision == 1);
    nTotal_cardinal   = nTotal_cardinal + numel(synthData_interleaved(idx).decision);

    nCorrect_oblique = 0; nTotal_oblique = 0;
    idx = strcmp({synthData_interleaved(:).image_task},'oblique') & cell2mat({synthData_interleaved(:).contrast_signed}) == 6;
    nCorrect_oblique = nCorrect_oblique + sum(synthData_interleaved(idx).decision == 2);
    nTotal_oblique   = nTotal_oblique + numel(synthData_interleaved(idx).decision);
    idx = strcmp({synthData_interleaved(:).image_task},'oblique') & cell2mat({synthData_interleaved(:).contrast_signed}) == -6;
    nCorrect_oblique = nCorrect_oblique + sum(synthData_interleaved(idx).decision == 1);
    nTotal_oblique   = nTotal_oblique + numel(synthData_interleaved(idx).decision);

    accuracy_cardinal(n) = nCorrect_cardinal / nTotal_cardinal;
    accuracy_oblique(n)  = nCorrect_oblique / nTotal_oblique;
    semAcc_cardinal(n)   = sqrt(accuracy_cardinal(n) .* (1 - accuracy_cardinal(n)) / nTotal_cardinal);
    semAcc_oblique(n)    = sqrt(accuracy_oblique(n) .* (1 - accuracy_oblique(n)) / nTotal_oblique);
end
%%
figure
set(gcf,'units','normalized','position',[0,0,1,1]);


for j1 = 1:nPrior
    for j2 = 1:nPrior
        subplot(nPrior, nPrior, (j1-1) * nPrior + j2); hold on
        idx_base = cardinal_prior == prior_list(j1) & oblique_prior == prior_list(j2);

        i_small = idx_base & cardinal_delta == delta_list(1) & oblique_delta == delta_list(1);
        errorbar(bPF(i_small),accuracy_cardinal(i_small),semAcc_cardinal(i_small),'LineWidth',2,'color','red','LineStyle','--');
        errorbar(bPF(i_small),accuracy_oblique(i_small),semAcc_oblique(i_small),'LineWidth',2,'color','blue','LineStyle','--');
        i_large = idx_base & cardinal_delta == delta_list(2) & oblique_delta == delta_list(2);
        errorbar(bPF(i_large),accuracy_cardinal(i_large),semAcc_cardinal(i_large),'LineWidth',2,'color','red','LineStyle','-');
        errorbar(bPF(i_large),accuracy_oblique(i_large),semAcc_oblique(i_large),'LineWidth',2,'color','blue','LineStyle','-');
        set(gca,'fontsize',16)
        box off
        xlabel('bPF')
        ylabel('Acc.')
        ylim([0.7,1])    
        
        title(sprintf('Cardinal prior = %.2f, Oblique prior = %.2f',prior_list(j1), prior_list(j2)))
    end
end
sgtitle('Accuracy at contrast 6','fontsize',20);
savename = fullfile(figFolder, 'Accuracy_contrast6_real_interleaved');
print(gcf, savename,'-dpng')
    


%% helper functions
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