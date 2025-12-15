clear all 
clc
close all

load('../../../results/eyeMovement/eyedata_summarized_interleaved');
%%%% This scripts visualize basic eye movement variables conditioned on
%%%% task context
save_folder = '../../../figures/figures_auto_interleavedpaper';
%%  Scatter plot of each variable, cardinal agains oblique task


eye_field_list = {'x_pos_avg';'y_pos_avg';'x_vel_avg';'y_vel_avg';'pupil_avg';...
                  'x_pos_var';'y_pos_var';'x_vel_var';'y_vel_var';'pupil_var'};


[t_values, p_values, df_values]  = deal(zeros(numel(eye_field_list), 2));
doPlot = 0;
for k = 1:numel(eye_field_list)
    field_name = eye_field_list{k};


    subjectCode = 'Ro';
    stats_info = plot_eye_scatter(eyeData_summarize_all, subjectCode, field_name, doPlot);
    t_values(k,1) = stats_info.tstat;
    p_values(k,1) = stats_info.p;
    df_values(k,1) = stats_info.df;

    subjectCode = 'Gr';
    stats_info = plot_eye_scatter(eyeData_summarize_all, subjectCode, field_name, doPlot);
    t_values(k,2) = stats_info.tstat;
    p_values(k,2) = stats_info.p;
    df_values(k,2) = stats_info.df;

end
rowNames    = eye_field_list;
colNames    = {'Monkey R';'Monkey G'};
save_name   = fullfile(save_folder,'eye_interleaved_compare_task.txt');
make_table(rowNames, colNames, t_values, p_values, df_values, save_name)

%%%%% only zero trials
eye_field_zero_list = {'x_pos_zero_avg';'y_pos_zero_avg';'x_vel_zero_avg';'y_vel_zero_avg';'pupil_zero_avg';...
                  'x_pos_zero_var';'y_pos_zero_var';'x_vel_zero_var';'y_vel_zero_var';'pupil_zero_var'};
[t_values_zero, p_values_zero, df_values_zero]  = deal(zeros(numel(eye_field_zero_list), 2));
for k = 1:numel(eye_field_zero_list)
    field_name = eye_field_zero_list{k};


    subjectCode = 'Ro';
    stats_info = plot_eye_scatter(eyeData_summarize_all, subjectCode, field_name, doPlot);
    t_values_zero(k,1) = stats_info.tstat;
    p_values_zero(k,1) = stats_info.p;
    df_values_zero(k,1) = stats_info.df;

    subjectCode = 'Gr';
    stats_info = plot_eye_scatter(eyeData_summarize_all, subjectCode, field_name, doPlot);
    t_values_zero(k,2) = stats_info.tstat;
    p_values_zero(k,2) = stats_info.p;
    df_values_zero(k,2) = stats_info.df;

end
rowNames    = eye_field_zero_list;
colNames    = {'Monkey R';'Monkey G'};
save_name   = fullfile(save_folder,'eye_interleaved_compare_task_zerotrials.txt');
make_table(rowNames, colNames, t_values_zero, p_values_zero, df_values_zero, save_name)
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
                  x_cardinal - x_cardinal_CI(:,1), x_cardinal_CI(:,2) - x_cardinal ,'.');
        box off
        xlabel('Cardinal trials');
        ylabel('Oblique trials');
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