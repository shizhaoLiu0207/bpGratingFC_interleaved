function h = plot_synth_interleaved_psycurve(synthData_interleaved, plotOptions)

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

xlabel('Contrast');
ylabel('Prob. Choose 0/45');
box off
line([0,0],[0,1],'color','black','linestyle','--');
line([-15,15],[0.5,0.5],'color','black','linestyle','--');
set(gca,'fontsize',plotOptions.ftsize)
legend(h, 'Cardinal task','Oblique task','location','southeast')
end