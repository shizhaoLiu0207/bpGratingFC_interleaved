clear all
clc
close all

global   bpGlobal  ftsize
bpGratingFCGlobal();
%%

figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,18,4])
flex_condition_list = {'Flexible';'Non-flexible'; 'Partially-flexible'};
nCondition = numel(flex_condition_list);
for n = 1:nCondition
    subplot(1,nCondition,n)
    flex_condition = flex_condition_list{n}; % switch, same
    
    x_pos = [0.5,5.5,7.5,2.5];
    
    switch flex_condition
        case 'Flexible'
            deltaI_cardinal = 1;
            deltaI_cardinal_cross = 0.2;
            
            deltaI_oblique = 1;
            deltaI_oblique_cross = 0.2;
        case 'Non-flexible'
            deltaI_cardinal = 0.6;
            deltaI_cardinal_cross = 0.6;
            
            deltaI_oblique = 0.6;
            deltaI_oblique_cross = 0.6;
        case 'Partially-flexible'
            deltaI_cardinal = 1;
            deltaI_cardinal_cross = 0.2;
            
            deltaI_oblique = 0.6;
            deltaI_oblique_cross = 0.6;
    end
    
    hold on
    
    bar(x_pos(1), deltaI_cardinal, 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_cardinal);
    bar(x_pos(2), deltaI_cardinal_cross, 'facecolor',bpGlobal.color_list.color_cardinal, 'EdgeColor',bpGlobal.color_list.color_oblique, 'linewidth',3);
    
    bar(x_pos(3), deltaI_oblique, 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_oblique);
    bar(x_pos(4), deltaI_oblique_cross, 'facecolor',bpGlobal.color_list.color_oblique, 'EdgeColor',bpGlobal.color_list.color_cardinal, 'linewidth',3);
    
    ylim([0,1.3])
    
    set(gca, 'fontsize', ftsize)
         set(gca, 'TickLabelInterpreter','tex')
    text(0.1,1.15,'Cardinal Context','color','red','FontSize',14,'FontWeight','bold')
    text(5.1,1.15,'Oblique Context','color','blue','FontSize',14,'FontWeight','bold')
    
    set(gca,'ytick',[],'xtick', [0.5,2.5,5.5,7.5], 'xticklabels', ...
                    {'\color{red}{Cardinal}';'\color{blue}{Oblique}'; '\color{red}{Cardinal}';'\color{blue}{Oblique}'});
    ylabel('$I_\textrm{redundancy}$','Interpreter','latex','FontSize',ftsize+4);
    title(flex_condition,'FontSize',ftsize+4);

end

save_folder = '../../figures/figures_final/fisher_info_bar';
save_name = fullfile(save_folder, sprintf('redundacy_illustration_bar.svg'));
print(save_name, '-dsvg');