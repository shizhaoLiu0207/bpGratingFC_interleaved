function plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
orientationsDEG = psyKernel_model(1).orientationsDEG;
 switch image_task
        case 'cardinal'
            plot_color = 'red';
            ideal_kernel            = 0.4*sin(pi * (orientationsDEG - 45)/90);
        case 'oblique'
            plot_color = 'blue';
            ideal_kernel            = 0.4*sin(pi * (orientationsDEG - 90)/90);
 end

switch kernel_type
    case 'spatial'
            hold on
            line([0,180],[0,0],'linestyle','--','color','black')
            plot(orientationsDEG,ideal_kernel,'linestyle','--','linewidth',2,'color',plot_color);
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

assert(numel(prior_task_list) == numel(mode_list), 'Number of prior and number of mode do not match')
for n = 1:numel(prior_task_list)
    prior_cardinal = prior_task_list{n}(1);
    prior_oblique = prior_task_list{n}(2);

    idx = find(cell2mat({psyKernel_model(:).prior_cardinal}) == prior_cardinal &...
          cell2mat({psyKernel_model(:).prior_oblique})  == prior_oblique & ...
          strcmp({psyKernel_model(:).image_task}, image_task) & ...
          strcmp({psyKernel_model(:).mode}, mode_list{n}));
    
   
    
    w_ori_all = {psyKernel_model(idx).w_ori};
    w_ori_all = cat(2,w_ori_all{:});
    
    w_time_all = {psyKernel_model(idx).w_time};
    w_time_all = cat(2,w_time_all{:});
    
    
    w_ori_avg = mean(w_ori_all,2);
    w_ori_sem = std(w_ori_all, [], 2) / sqrt(size(w_ori_all, 2));
    
    w_time_avg = mean(w_time_all,2);
    w_time_sem = std(w_time_all, [], 2) / sqrt(size(w_time_all, 2));
    
   
    switch kernel_type
        case 'spatial'
            plot(orientationsDEG, w_ori_avg,'Color',plot_color,'LineWidth',2); hold on
            fill([orientationsDEG, fliplr(orientationsDEG)], [w_ori_avg' - w_ori_sem', fliplr(w_ori_avg' + w_ori_sem')],...
                plot_color, 'FaceAlpha',0.5);
            
        case 'temporal'
            x = [1:numel(w_time_avg)];
            
            plot(x, w_time_avg,'Color',plot_color,'LineWidth',2);  hold on
            fill([x, fliplr(x)], [w_time_avg' - w_time_sem', fliplr(w_time_avg' + w_time_sem')],...
                plot_color, 'FaceAlpha',0.5);
    end
    set(gca,'fontsize',18);
    box off
end
end