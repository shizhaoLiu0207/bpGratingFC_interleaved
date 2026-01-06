function [r_cardinal, r_oblique, nSession] = plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode)
orientationsDEG = psyKernel_model(1).orientationsDEG;
ideal_kernel_cardinal            = 0.4*sin(pi * (orientationsDEG - 45)/90);
ideal_kernel_oblique             = 0.4*sin(pi * (orientationsDEG - 90)/90);
 switch image_task
        case 'cardinal'
            plot_color = 'red';
            ideal_kernel_plot = ideal_kernel_cardinal;
        case 'oblique'
            plot_color = 'blue';
            ideal_kernel_plot = ideal_kernel_oblique;
            
 end

switch kernel_type
    case 'spatial'
            hold on
            line([0,180],[0,0],'linestyle','--','color','black')
            plot(orientationsDEG,ideal_kernel_plot,'linestyle','--','linewidth',2,'color',plot_color);
            line([0,0],[-0.6,0.6],'linestyle','--','color','black'); 
            line([90,90],[-0.6,0.6],'linestyle','--','color','black');
            line([45,45],[-0.6,0.6],'linestyle','--','color','black');
            line([135,135],[-0.6,0.6],'linestyle','--','color','black');
            xlim([-10,190])
            ylim([-0.5,0.5])
            xlabel('Orientation');
            ylabel('Spatial kernel')
    case 'temporal'
        hold on
end



prior_cardinal = prior_task(1);
prior_oblique = prior_task(2);

idx = find(cell2mat({psyKernel_model(:).prior_cardinal}) == prior_cardinal &...
      cell2mat({psyKernel_model(:).prior_oblique})  == prior_oblique & ...
      strcmp({psyKernel_model(:).image_task}, image_task) & ...
      strcmp({psyKernel_model(:).mode}, mode));



w_ori_all = {psyKernel_model(idx).w_ori};
w_ori_all = cat(2,w_ori_all{:});

w_time_all = {psyKernel_model(idx).w_time};
w_time_all = cat(2,w_time_all{:});


w_ori_avg = mean(w_ori_all,2);
w_ori_sem = std(w_ori_all, [], 2) / sqrt(size(w_ori_all, 2));

w_time_avg = mean(w_time_all,2);
w_time_sem = std(w_time_all, [], 2) / sqrt(size(w_time_all, 2));

%%%% correlation between the two ideal kernels
r_cardinal  = corr(w_ori_avg, ideal_kernel_cardinal');
r_oblique   = corr(w_ori_avg, ideal_kernel_oblique');  

nSession    = numel(idx); 

switch kernel_type
    case 'spatial'
        plot(orientationsDEG, w_ori_avg,'Color',plot_color,'LineWidth',2); hold on
        fill([orientationsDEG, fliplr(orientationsDEG)], [w_ori_avg' - w_ori_sem', fliplr(w_ori_avg' + w_ori_sem')],...
            plot_color, 'FaceAlpha',0.5);
        % for t = 1:size(w_ori_all,2)
        %     plot(orientationsDEG, w_ori_all(:,t),'Color',plot_color,'LineWidth',0.5); 
        % end
    case 'temporal'
        x = [1:numel(w_time_avg)];
        
        plot(x, w_time_avg,'Color',plot_color,'LineWidth',2);  hold on
        fill([x, fliplr(x)], [w_time_avg' - w_time_sem', fliplr(w_time_avg' + w_time_sem')],...
            plot_color, 'FaceAlpha',0.5);
        % for t = 1:size(w_ori_all,2)
        %     plot(orientationsDEG, w_ori_all(:,t),'Color',plot_color,'LineWidth',0.5); 
        % end
end
set(gca,'fontsize',18);
box off

end