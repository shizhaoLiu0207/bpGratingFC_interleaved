clear all 
clc
close all

load('../../../results/eyeMovement/eyedata_summarized_interleaved');
%%%% This scripts visualize basic eye movement variables conditioned on
%%%% task context
save_folder = '../../../figures/figures_auto_interleavedpaper';


eyeData_summarize_all = normalize_pupil(eyeData_summarize_all);

idx_rm = ismember({eyeData_summarize_all(:).sessionStr},{'Ro20220107';'Ro20220220'});
eyeData_summarize_all(idx_rm) = [];

%%  Scatter plot of each variable, cardinal agains oblique task

use_zero_only = true;
if use_zero_only
    eye_field_list = {'x_pos_zero_avg';'y_pos_zero_avg';'x_vel_zero_avg';'y_vel_zero_avg';'pupil_zero_avg';...
                      'x_pos_zero_var';'y_pos_zero_var';'x_vel_zero_var';'y_vel_zero_var';'pupil_zero_var'};
    save_name   = fullfile(save_folder,'eye_interleaved_compare_task_zerotrials.svg');
else
    eye_field_list = {'x_pos_avg';'y_pos_avg';'x_vel_avg';'y_vel_avg';'pupil_avg';...
                      'x_pos_var';'y_pos_var';'x_vel_var';'y_vel_var';'pupil_var'};

    save_name   = fullfile(save_folder,'eye_interleaved_compare_task.svg');
end
subplot_index_list = [1,5,9,13,17,2,6,10,14,18];

[t_values, p_values, df_values]  = deal(zeros(numel(eye_field_list), 2));
doPlot = 1;
figure;
set(gcf, 'unit','normalized','position',[0,0,1,1]);

for k = 1:numel(eye_field_list)
    field_name = eye_field_list{k};
   

    subjectCode = 'Ro';
  
    subplot(5,4, subplot_index_list(k));

    stats_info = plot_eye_scatter(eyeData_summarize_all, subjectCode, field_name, doPlot);
    if stats_info.p < 0.05
        title(sprintf('\\textbf{%s}, \\boldmath $t(%d) = %.2f^{%s}$',strrep(field_name,'_','-'), stats_info.df, stats_info.tstat, fig.p2star(stats_info.p))...
            ,'Interpreter','latex','FontSize',16);  
    else
        title(sprintf('%s, $t(%d) = %.2f^{%s}$',strrep(field_name,'_','-'), stats_info.df, stats_info.tstat, fig.p2star(stats_info.p))...
            ,'Interpreter','latex','FontSize',16);  
    end
    % t_values(k,1) = stats_info.tstat;
    % p_values(k,1) = stats_info.p;
    % df_values(k,1) = stats_info.df;
  

    subjectCode = 'Gr';
    
    subplot(5,4, subplot_index_list(k) + 2);
    stats_info = plot_eye_scatter(eyeData_summarize_all, subjectCode, field_name, doPlot);

    if stats_info.p < 0.05
        title(sprintf('\\textbf{%s}, \\boldmath $t(%d) = %.2f^{%s}$',strrep(field_name,'_','-'), stats_info.df, stats_info.tstat, fig.p2star(stats_info.p))...
            ,'Interpreter','latex','FontSize',16);  
    else
        title(sprintf('%s, $t(%d) = %.2f^{%s}$',strrep(field_name,'_','-'), stats_info.df, stats_info.tstat, fig.p2star(stats_info.p))...
            ,'Interpreter','latex','FontSize',16);  
    end


   
end



annotation('textbox',[0.28,0.95,0.15,0.04],'string','Monkey R','fontsize',20,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.69,0.95,0.15,0.04],'string','Monkey G','fontsize',20,'FontWeight','bold','EdgeColor','none')


print(save_name,'-dsvg','-vector');

%%
function stats_info = plot_eye_scatter(eyeData_summarize_all, subjectCode, field_name, doPlot)

    idx_session = contains({eyeData_summarize_all(:).sessionStr},subjectCode);
    
    eval(sprintf('x       = {eyeData_summarize_all(idx_session).%s_median};', field_name));
    eval(sprintf('x_CI    = {eyeData_summarize_all(idx_session).%s_CI};', field_name)); 
    
    x       = cat(1,x{:}); 
    x_CI    = cat(3,x_CI{:}); 
    x_cardinal  = x(:,1);
    x_oblique   = x(:,2);
    x_cardinal_CI = squeeze(x_CI(:,1,:))';
    x_oblique_CI  = squeeze(x_CI(:,2,:))';
    if doPlot
        errorbar(x_cardinal, x_oblique, ...
                  x_oblique - x_oblique_CI(:,1), x_oblique_CI(:,2) -  x_oblique,...
                  x_cardinal - x_cardinal_CI(:,1), x_cardinal_CI(:,2) - x_cardinal ,'.','LineWidth',2);
        box off
        xlabel('Cardinal trials','color','red','FontSize',14);
        ylabel('Oblique trials','color','blue','FontSize',14);
        xmin = min(x_cardinal_CI(:));
        xmax = max(x_cardinal_CI(:));
        ymin = min(x_oblique_CI(:));
        ymax = max(x_oblique_CI(:));
        line([xmin,xmax], [ymin, ymax], 'linewidth',1.5,'color','black','linestyle','--')
    end
    [~,p,~,stats] = ttest(x_cardinal, x_oblique);
    stats_info.tstat = stats.tstat;
    stats_info.p = p;
    stats_info.df = stats.df;
end






function make_table(rowNames, colNames, t_values, p_values, df_values, save_name)

    %rowNames = eye_field_list;
    for k = 1:numel(rowNames)
        rowNames{k} = replace(rowNames{k},'_','-'); 
    end
    %colNames = {'MonkeyR, cardinal';'MonkeyR, oblique';'MonkeyG, oblique';'MonkeyG, cardinal'};
    
    formatted_cells = strings(size(t_values));
    for i = 1:size(t_values, 1)
        for j = 1:size(t_values, 2)
            if p_values(i,j) < 0.05
                formatted_cells(i,j) = sprintf('\\makecell{$\\mathbf{t(%d)=%.2f}$, \\\\ $p=\\num{%.2e}$}', df_values(i,j), t_values(i,j), p_values(i,j));
            else
                formatted_cells(i,j) = sprintf('\\makecell{$t(%d)=%.2f$, \\\\ $p=\\num{%.2e}$}', df_values(i,j), t_values(i,j), p_values(i,j));
            end
        end
    end
    
    % Convert to table
    T = array2table(formatted_cells, 'VariableNames', colNames, 'RowNames', rowNames);
    
    % Export to LaTeX-friendly format
    % We'll write to file using custom delimiters and without quotes
    %filename = fullfile(save_folder,'figureS21_eyeMetric_partial_correlation_tabel.txt');
    fid = fopen(save_name, 'w');
    fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('c', 1, width(T)));
    
    % Header
    fprintf(fid, ' & %s \\\\\n\\midrule\n', strjoin(T.Properties.VariableNames, ' & '));
    
    % Rows
    for i = 1:height(T)
        rowLabel = T.Properties.RowNames{i};
        rowData = T{i,:};
        fprintf(fid, '%s & %s \\\\\n', rowLabel, strjoin(rowData, ' & '));
           if i < height(T)
           fprintf(fid, '\\midrule\n');  % line between rows
            
        end
    end
    
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
    fclose(fid);

end

function eyeData_summarize_all = normalize_pupil(eyeData_summarize_all)

    idx = contains({eyeData_summarize_all(:).sessionStr}, 'Ro');
    raw_pupil =[eyeData_summarize_all(idx).pupil_avg_median];
    raw_pupil_mean_rolo = mean(raw_pupil(:));
    raw_pupil_std_rolo  = std(raw_pupil(:));
    
    idx = contains({eyeData_summarize_all(:).sessionStr}, 'Gr');
    raw_pupil = [eyeData_summarize_all(idx).pupil_avg_median];
    raw_pupil_mean_gremlin = mean(raw_pupil(:)); 
    raw_pupil_std_gremlin  = std(raw_pupil(:));
    
    for i = 1:numel(eyeData_summarize_all)
        switch eyeData_summarize_all(i).sessionStr(1:2)
            case 'Ro'
                eyeData_summarize_all(i).pupil_avg_median = ...
                    (eyeData_summarize_all(i).pupil_avg_median - raw_pupil_mean_rolo) / raw_pupil_std_rolo;
                eyeData_summarize_all(i).pupil_avg_CI = ...
                     (eyeData_summarize_all(i).pupil_avg_CI - raw_pupil_mean_rolo) / raw_pupil_std_rolo;
                eyeData_summarize_all(i).pupil_var_median = ...
                    eyeData_summarize_all(i).pupil_var_median / raw_pupil_std_rolo^2;
                eyeData_summarize_all(i).pupil_var_CI = ...
                    eyeData_summarize_all(i).pupil_var_CI / raw_pupil_std_rolo^2;

                eyeData_summarize_all(i).pupil_zero_avg_median = ...
                    (eyeData_summarize_all(i).pupil_zero_avg_median - raw_pupil_mean_rolo) / raw_pupil_std_rolo;
                eyeData_summarize_all(i).pupil_zero_avg_CI = ...
                     (eyeData_summarize_all(i).pupil_zero_avg_CI - raw_pupil_mean_rolo) / raw_pupil_std_rolo;
                eyeData_summarize_all(i).pupil_zero_var_median = ...
                    eyeData_summarize_all(i).pupil_zero_var_median / raw_pupil_std_rolo^2;
                eyeData_summarize_all(i).pupil_zero_var_CI = ...
                    eyeData_summarize_all(i).pupil_zero_var_CI / raw_pupil_std_rolo^2;
            case 'Gr'
                eyeData_summarize_all(i).pupil_avg_median = ...
                    (eyeData_summarize_all(i).pupil_avg_median - raw_pupil_mean_gremlin) / raw_pupil_std_gremlin;
                eyeData_summarize_all(i).pupil_avg_CI = ...
                     (eyeData_summarize_all(i).pupil_avg_CI - raw_pupil_mean_gremlin) / raw_pupil_std_gremlin;
                eyeData_summarize_all(i).pupil_var_median = ...
                    eyeData_summarize_all(i).pupil_var_median / raw_pupil_std_gremlin^2;
                eyeData_summarize_all(i).pupil_var_CI = ...
                    eyeData_summarize_all(i).pupil_var_CI / raw_pupil_std_gremlin^2;

                eyeData_summarize_all(i).pupil_zero_avg_median = ...
                    (eyeData_summarize_all(i).pupil_zero_avg_median - raw_pupil_mean_gremlin) / raw_pupil_std_gremlin;
                eyeData_summarize_all(i).pupil_zero_avg_CI = ...
                     (eyeData_summarize_all(i).pupil_zero_avg_CI - raw_pupil_mean_gremlin) / raw_pupil_std_gremlin;
                eyeData_summarize_all(i).pupil_zero_var_median = ...
                    eyeData_summarize_all(i).pupil_zero_var_median / raw_pupil_std_gremlin^2;
                eyeData_summarize_all(i).pupil_zero_var_CI = ...
                    eyeData_summarize_all(i).pupil_zero_var_CI / raw_pupil_std_gremlin^2;

        end
    end


end