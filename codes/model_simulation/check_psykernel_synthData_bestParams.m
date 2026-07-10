clear all
clc
close all


global bpGlobal
bpGratingFCGlobal()
orientationsDEG                  = [0:15:165];
w_ori_cardinal_ideal            = 0.4*sin(pi * (orientationsDEG - 45)/90)';
w_ori_oblique_ideal             = 0.4*sin(pi * (orientationsDEG - 90)/90)';

%% 
w_ori_cardinal_ideal_multiple = 0.5* repmat(w_ori_cardinal_ideal, 1, 50);
w_ori_oblique_ideal_multiple = repmat(w_ori_oblique_ideal, 1, 50);

w_ori_cardinal_ideal_multiple = w_ori_cardinal_ideal_multiple + 0.2 * randn(size(w_ori_cardinal_ideal_multiple));
w_ori_oblique_ideal_multiple  = w_ori_oblique_ideal_multiple + 0.2 * randn(size(w_ori_oblique_ideal_multiple));

[dist_avg, dist_combination] = compute_distance_kernels(w_ori_cardinal_ideal_multiple, w_ori_oblique_ideal_multiple, w_ori_cardinal_ideal);
%% check results
%%%%%% Compute correlations between kernels
%%%%% load empirical kernel

psyKernel_Ro = load('/Users/liushizhao/projects_local/bpGratingEx/results/behavior/Rolo_psyKernel_table_final');
psyKernel_Gr = load('/Users/liushizhao/projects_local/bpGratingEx/results/behavior/Gremlin_psyKernel_table');

sessionName_list  = bpGlobal.rolo.session_list.switching;
idx_list          = ismember({psyKernel_Ro.psyKernel_table(:).sessionName},sessionName_list);

tmp                         = {psyKernel_Ro.psyKernel_table(idx_list).w_ori_cardinal};
w_ori_cardinal_Ro = cat(2,tmp{:});
tmp                         = {psyKernel_Ro.psyKernel_table(idx_list).w_ori_oblique};
w_ori_oblique_Ro  = cat(2,tmp{:});

sessionName_list  = bpGlobal.gremlin.session_list.interleaved_real;
idx_list          = ismember({psyKernel_Gr.psyKernel_table(:).sessionName},sessionName_list);

tmp                         = {psyKernel_Gr.psyKernel_table(idx_list).w_ori_cardinal};
w_ori_cardinal_Gr = cat(2,tmp{:});
tmp                         = {psyKernel_Gr.psyKernel_table(idx_list).w_ori_oblique};
w_ori_oblique_Gr  = cat(2,tmp{:});





%compute_assymetry_kernels(w_ori_cardinal_ideal/2, w_ori_oblique_ideal, w_ori_cardinal_ideal);

subject_list = {'monkeyR';'monkeyG'};
data_folder  = '../../results/model_behavior/psyKernel_bestParams';
psyKernel_corr = struct();
%psyKernel_dist = struct();
k  = 1;
for i = 1:2
    for i_set = 1:10
        subjectCode = subject_list{i};
   
        data_name = fullfile(data_folder, sprintf('psyKernels_best_parameters_%s_set_%d', subjectCode, i_set));
        load(data_name);
        
        
        idx_cardinal = strcmp({psyKernel_model(:).task},'cardinal');
        idx_oblique = strcmp({psyKernel_model(:).task},'oblique');
        prior_cardinal  =  psyKernel_model(find(idx_cardinal,1)).prior;
        prior_oblique   =  psyKernel_model(find(idx_oblique,1)).prior;
        delta_cardinal  =  psyKernel_model(find(idx_cardinal,1)).delta;
        delta_oblique   = psyKernel_model(find(idx_oblique,1)).delta;
        tmp = {psyKernel_model(idx_cardinal).w_ori};
        w_ori_cardinal_model = cat(2, tmp{:});
        tmp = {psyKernel_model(idx_oblique).w_ori};
        w_ori_oblique_model = cat(2, tmp{:});
        
        switch subjectCode
            case 'monkeyR'
                w_ori_cardinal_empirical = w_ori_cardinal_Ro;
                w_ori_oblique_empirical = w_ori_oblique_Ro;
            case 'monkeyG'
                w_ori_cardinal_empirical = w_ori_cardinal_Gr;
                w_ori_oblique_empirical = w_ori_oblique_Gr;
        end
        %%% corr_model_empirical
        [r_model_empirical_cardinal_avg, r_model_empirical_cardinal_combinations] = ...
                compute_corr_kernels(w_ori_cardinal_model, w_ori_cardinal_empirical);
        [r_model_empirical_oblique_avg, r_model_empirical_oblique_combinations] = ...
                compute_corr_kernels(w_ori_oblique_model, w_ori_oblique_empirical);
        
        %%% corr_model_ideal
        [r_model_ideal_cardinal_avg, r_model_ideal_cardinal_combinations] = ...
                compute_corr_kernels(w_ori_cardinal_model, w_ori_cardinal_ideal);
        [r_model_ideal_oblique_avg, r_model_ideal_oblique_combinations] = ...
                compute_corr_kernels(w_ori_oblique_model, w_ori_oblique_ideal);

        %%% corr_empirical_ideal
        [r_empirical_ideal_cardinal_avg, ~] = ...
                compute_corr_kernels(w_ori_cardinal_empirical, w_ori_cardinal_ideal);
        [r_empirical_ideal_oblique_avg, ~] = ...
                compute_corr_kernels(w_ori_oblique_empirical, w_ori_oblique_ideal);


        
        psyKernel_corr(k).subjectCode                               = subjectCode;
        psyKernel_corr(k).i_set                                     = i_set;
        psyKernel_corr(k).prior_cardinal                            = prior_cardinal;
        psyKernel_corr(k).prior_oblique                             = prior_oblique;
        psyKernel_corr(k).delta_cardinal                            = delta_cardinal;
        psyKernel_corr(k).delta_oblique                             = delta_oblique;
        psyKernel_corr(k).asymmetry_prior               = prior_cardinal - prior_oblique;
        psyKernel_corr(k).asymmetry_delta               = delta_cardinal - delta_oblique;
        psyKernel_corr(k).w_cardinal_model                          = mean(w_ori_cardinal_model, 2);
        psyKernel_corr(k).w_oblique_model                           = mean(w_ori_oblique_model, 2);
        psyKernel_corr(k).w_cardinal_empirical                      = mean(w_ori_cardinal_empirical, 2);
        psyKernel_corr(k).w_oblique_empirical                       = mean(w_ori_oblique_empirical, 2);
        psyKernel_corr(k).w_cardinal_ideal                           = w_ori_cardinal_ideal;
        psyKernel_corr(k).w_oblique_ideal                            = w_ori_oblique_ideal;
        
        psyKernel_corr(k).r_model_empirical_cardinal_avg            = r_model_empirical_cardinal_avg;
        psyKernel_corr(k).r_model_empirical_cardinal_combinations   = r_model_empirical_cardinal_combinations;
        
        psyKernel_corr(k).r_model_empirical_oblique_avg             = r_model_empirical_oblique_avg;
        psyKernel_corr(k).r_model_empirical_oblique_combinations    = r_model_empirical_oblique_combinations;
        
        psyKernel_corr(k).r_model_ideal_cardinal_avg                = r_model_ideal_cardinal_avg;
        psyKernel_corr(k).r_model_ideal_cardinal_combinations       = r_model_ideal_cardinal_combinations;
        psyKernel_corr(k).r_model_ideal_oblique_avg                 = r_model_ideal_oblique_avg;
        psyKernel_corr(k).r_model_ideal_oblique_combinations        = r_model_ideal_oblique_combinations;

        psyKernel_corr(k).r_empirical_ideal_cardinal_avg            = r_empirical_ideal_cardinal_avg;
        psyKernel_corr(k).r_empirical_ideal_oblique_avg             = r_empirical_ideal_oblique_avg;


        
        
         %%% distance_model_empirical
        [dist_model_empirical_cardinal_avg, dist_model_empirical_cardinal_combinations] = ...
                compute_distance_kernels(w_ori_cardinal_model, w_ori_cardinal_empirical, w_ori_cardinal_ideal);
        [dist_model_empirical_oblique_avg, dist_model_empirical_oblique_combinations] = ...
                compute_distance_kernels(w_ori_oblique_model, w_ori_oblique_empirical, w_ori_cardinal_ideal);
        
        %%% distance_model_ideal
        [dist_model_ideal_cardinal_avg, dist_model_ideal_cardinal_combinations] = ...
                compute_distance_kernels(w_ori_cardinal_model, w_ori_cardinal_ideal, w_ori_cardinal_ideal);
        [dist_model_ideal_oblique_avg, dist_model_ideal_oblique_combinations] = ...
                compute_distance_kernels(w_ori_oblique_model, w_ori_oblique_ideal, w_ori_cardinal_ideal);

        %%% distance_empirical_ideal
        [dist_empirical_ideal_cardinal_avg, dist_empirical_ideal_cardinal_combinations] = ...
                compute_distance_kernels(w_ori_cardinal_empirical, w_ori_cardinal_ideal, w_ori_cardinal_ideal);
        [dist_empirical_ideal_oblique_avg, dist_empirical_ideal_oblique_combinations] = ...
                compute_distance_kernels(w_ori_oblique_empirical, w_ori_oblique_ideal, w_ori_cardinal_ideal);



        %%%% a distance-based metric for assymetry
        [asymmetry_model_avg, asymmetry_model_combinations] =...
            compute_distance_kernels(w_ori_cardinal_model, w_ori_oblique_model, w_ori_cardinal_ideal);
        
       

        psyKernel_corr(k).dist_model_empirical_cardinal_avg = dist_model_empirical_cardinal_avg;
        psyKernel_corr(k).dist_model_empirical_cardinal_combinations = dist_model_empirical_cardinal_combinations;
        psyKernel_corr(k).dist_model_empirical_oblique_avg = dist_model_empirical_oblique_avg;
        psyKernel_corr(k).dist_model_empirical_oblique_combinations = dist_model_empirical_oblique_combinations;
        psyKernel_corr(k).dist_model_ideal_cardinal_avg = dist_model_ideal_cardinal_avg;
        psyKernel_corr(k).dist_model_ideal_cardinal_combinations = dist_model_ideal_cardinal_combinations;
        psyKernel_corr(k).dist_model_ideal_oblique_avg = dist_model_ideal_oblique_avg;
        psyKernel_corr(k).dist_model_ideal_oblique_combinations = dist_model_ideal_oblique_combinations;
        psyKernel_corr(k).dist_empirical_ideal_cardinal_avg = dist_empirical_ideal_cardinal_avg;
        psyKernel_corr(k).dist_empirical_ideal_cardinal_combinations = dist_empirical_ideal_cardinal_combinations;
        psyKernel_corr(k).dist_empirical_ideal_oblique_avg = dist_empirical_ideal_oblique_avg;
        psyKernel_corr(k).dist_empirical_ideal_oblique_combinations = dist_empirical_ideal_oblique_combinations;


        psyKernel_corr(k).asymmetry_model_avg                            = asymmetry_model_avg;
        psyKernel_corr(k).asymmetry_model_combinations                   = asymmetry_model_combinations;
        

        k = k+1;
    end
end
%%% load results of psykernel of symmatric parameters



data_folder  = '../../results/model_behavior/psyKernel_baselineParams';
subject_list = {'monkeyR';'monkeyG'};
psyKernel_baseline_corr = struct();
%psyKernel_baseline_dist = struct();
k = 1;
for i = 1:2
    subjectCode = subject_list{i};
    for i_set = 1:30
     
    
        data_name = fullfile(data_folder, sprintf('psyKernels_baseline_parameters_set_%d', i_set));
        load(data_name);
        
        
        idx_cardinal = strcmp({psyKernel_model(:).task},'cardinal');
        idx_oblique = strcmp({psyKernel_model(:).task},'oblique');
        prior_cardinal  =  psyKernel_model(find(idx_cardinal,1)).prior;
        prior_oblique   =  psyKernel_model(find(idx_oblique,1)).prior;
        delta_cardinal  =  psyKernel_model(find(idx_cardinal,1)).delta;
        delta_oblique   = psyKernel_model(find(idx_oblique,1)).delta;
        tmp = {psyKernel_model(idx_cardinal).w_ori};
        w_ori_cardinal_model = cat(2, tmp{:});
        tmp = {psyKernel_model(idx_oblique).w_ori};
        w_ori_oblique_model = cat(2, tmp{:});
        
        switch subjectCode
            case 'monkeyR'
                w_ori_cardinal_empirical = w_ori_cardinal_Ro;
                w_ori_oblique_empirical = w_ori_oblique_Ro;
            case 'monkeyG'
                w_ori_cardinal_empirical = w_ori_cardinal_Gr;
                w_ori_oblique_empirical = w_ori_oblique_Gr;
        end
        %%% corr_model_empirical
        [r_model_empirical_cardinal_avg, r_model_empirical_cardinal_combinations] = ...
                compute_corr_kernels(w_ori_cardinal_model, w_ori_cardinal_empirical);
        [r_model_empirical_oblique_avg, r_model_empirical_oblique_combinations] = ...
                compute_corr_kernels(w_ori_oblique_model, w_ori_oblique_empirical);
        
        %%% corr_model_ideal
        [r_model_ideal_cardinal_avg, r_model_ideal_cardinal_combinations] = ...
                compute_corr_kernels(w_ori_cardinal_model, w_ori_cardinal_ideal);
        [r_model_ideal_oblique_avg, r_model_ideal_oblique_combinations] = ...
                compute_corr_kernels(w_ori_oblique_model, w_ori_oblique_ideal);
    
        %%% corr_empirical_ideal
        [r_empirical_ideal_cardinal_avg, ~] = ...
                compute_corr_kernels(w_ori_cardinal_empirical, w_ori_cardinal_ideal);
        [r_empirical_ideal_oblique_avg, ~] = ...
                compute_corr_kernels(w_ori_oblique_empirical, w_ori_oblique_ideal);
    

        
        psyKernel_baseline_corr(k).subjectCode                               = subjectCode;
        psyKernel_baseline_corr(k).i_set                                     = i_set;
        psyKernel_baseline_corr(k).prior_cardinal                            = prior_cardinal;
        psyKernel_baseline_corr(k).prior_oblique                             = prior_oblique;
        psyKernel_baseline_corr(k).delta_cardinal                            = delta_cardinal;
        psyKernel_baseline_corr(k).delta_oblique                             = delta_oblique;
        psyKernel_baseline_corr(k).asymmetry_prior               = prior_cardinal - prior_oblique;
        psyKernel_baseline_corr(k).asymmetry_delta               = delta_cardinal - delta_oblique;
        psyKernel_baseline_corr(k).w_cardinal_model                          = mean(w_ori_cardinal_model, 2);
        psyKernel_baseline_corr(k).w_oblique_model                           = mean(w_ori_oblique_model, 2);
        psyKernel_baseline_corr(k).w_cardinal_empirical                      = mean(w_ori_cardinal_empirical, 2);
        psyKernel_baseline_corr(k).w_oblique_empirical                       = mean(w_ori_oblique_empirical, 2);
        psyKernel_baseline_corr(k).w_cardinal_ideal                           = w_ori_cardinal_ideal;
        psyKernel_baseline_corr(k).w_oblique_ideal                            = w_ori_oblique_ideal;
        psyKernel_baseline_corr(k).r_model_empirical_cardinal_avg            = r_model_empirical_cardinal_avg;
        psyKernel_baseline_corr(k).r_model_empirical_cardinal_combinations   = r_model_empirical_cardinal_combinations;
        
        psyKernel_baseline_corr(k).r_model_empirical_oblique_avg             = r_model_empirical_oblique_avg;
        psyKernel_baseline_corr(k).r_model_empirical_oblique_combinations    = r_model_empirical_oblique_combinations;
        
        psyKernel_baseline_corr(k).r_model_ideal_cardinal_avg                = r_model_ideal_cardinal_avg;
        psyKernel_baseline_corr(k).r_model_ideal_cardinal_combinations       = r_model_ideal_cardinal_combinations;
        psyKernel_baseline_corr(k).r_model_ideal_oblique_avg                 = r_model_ideal_oblique_avg;
        psyKernel_baseline_corr(k).r_model_ideal_oblique_combinations        = r_model_ideal_oblique_combinations;
    
        psyKernel_baseline_corr(k).r_empirical_ideal_cardinal_avg            = r_empirical_ideal_cardinal_avg;
        psyKernel_baseline_corr(k).r_empirical_ideal_oblique_avg             = r_empirical_ideal_oblique_avg;
    


         %%% distance_model_empirical
        [dist_model_empirical_cardinal_avg, dist_model_empirical_cardinal_combinations] = ...
                compute_distance_kernels(w_ori_cardinal_model, w_ori_cardinal_empirical, w_ori_cardinal_ideal);
        [dist_model_empirical_oblique_avg, dist_model_empirical_oblique_combinations] = ...
                compute_distance_kernels(w_ori_oblique_model, w_ori_oblique_empirical, w_ori_cardinal_ideal);
        
        %%% distance_model_ideal
        [dist_model_ideal_cardinal_avg, dist_model_ideal_cardinal_combinations] = ...
                compute_distance_kernels(w_ori_cardinal_model, w_ori_cardinal_ideal, w_ori_cardinal_ideal);
        [dist_model_ideal_oblique_avg, dist_model_ideal_oblique_combinations] = ...
                compute_distance_kernels(w_ori_oblique_model, w_ori_oblique_ideal, w_ori_cardinal_ideal);

        %%% distance_empirical_ideal
        [dist_empirical_ideal_cardinal_avg, dist_empirical_ideal_cardinal_combinations] = ...
                compute_distance_kernels(w_ori_cardinal_empirical, w_ori_cardinal_ideal, w_ori_cardinal_ideal);
        [dist_empirical_ideal_oblique_avg, dist_empirical_ideal_oblique_combinations] = ...
                compute_distance_kernels(w_ori_oblique_empirical, w_ori_oblique_ideal, w_ori_cardinal_ideal);

       %%%% a distance-based metric for assymetry
        [asymmetry_model_avg, asymmetry_model_combinations] =...
            compute_distance_kernels(w_ori_cardinal_model, w_ori_oblique_model, w_ori_cardinal_ideal);
        
       
        psyKernel_baseline_corr(k).dist_model_empirical_cardinal_avg = dist_model_empirical_cardinal_avg;
        psyKernel_baseline_corr(k).dist_model_empirical_cardinal_combinations = dist_model_empirical_cardinal_combinations;
        psyKernel_baseline_corr(k).dist_model_empirical_oblique_avg = dist_model_empirical_oblique_avg;
        psyKernel_baseline_corr(k).dist_model_empirical_oblique_combinations = dist_model_empirical_oblique_combinations;
        psyKernel_baseline_corr(k).dist_model_ideal_cardinal_avg = dist_model_ideal_cardinal_avg;
        psyKernel_baseline_corr(k).dist_model_ideal_cardinal_combinations = dist_model_ideal_cardinal_combinations;
        psyKernel_baseline_corr(k).dist_model_ideal_oblique_avg = dist_model_ideal_oblique_avg;
        psyKernel_baseline_corr(k).dist_model_ideal_oblique_combinations = dist_model_ideal_oblique_combinations;
        psyKernel_baseline_corr(k).dist_empirical_ideal_cardinal_avg = dist_empirical_ideal_cardinal_avg;
        psyKernel_baseline_corr(k).dist_empirical_ideal_cardinal_combinations = dist_empirical_ideal_cardinal_combinations;
        psyKernel_baseline_corr(k).dist_empirical_ideal_oblique_avg = dist_empirical_ideal_oblique_avg;
        psyKernel_baseline_corr(k).dist_empirical_ideal_oblique_combinations = dist_empirical_ideal_oblique_combinations;

        psyKernel_baseline_corr(k).asymmetry_model_avg                            = asymmetry_model_avg;
        psyKernel_baseline_corr(k).asymmetry_model_combinations                   = asymmetry_model_combinations;

    
        k = k+1;
    end
end


%% visualization: model-empirical
plotOptions = 'corr'; % corr or dist
doErrorbar = true;
%%%%% scatter, corr_model_empirical_cardinal, corr_model_empirical_oblique
%%%%% scatter, corr_model_ideal_cardinal, corr_model_ideal_oblique. With
%%%%% the corr_empirical_ideal marked
colors_list = get(groot, 'defaultAxesColorOrder');
figure
for i  = 1:2
    subjectCode = subject_list{i};

    subplot(1,2,i); hold on

    idx = strcmp({psyKernel_corr(:).subjectCode}, subjectCode);
    switch plotOptions
        case 'corr'
            x_avg = [psyKernel_corr(idx).r_model_empirical_cardinal_avg];
            y_avg = [psyKernel_corr(idx).r_model_empirical_oblique_avg];
            x_sem = cellfun(@std,{psyKernel_corr(idx).r_model_empirical_cardinal_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).r_model_empirical_cardinal_combinations}));
            y_sem = cellfun(@std,{psyKernel_corr(idx).r_model_empirical_oblique_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).r_model_empirical_oblique_combinations}));
        case 'dist'
            x_avg = [psyKernel_corr(idx).dist_model_empirical_cardinal_avg];
            y_avg = [psyKernel_corr(idx).dist_model_empirical_oblique_avg];
            x_sem = cellfun(@std,{psyKernel_corr(idx).dist_model_empirical_cardinal_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).dist_model_empirical_cardinal_combinations}));
            y_sem = cellfun(@std,{psyKernel_corr(idx).dist_model_empirical_oblique_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).dist_model_empirical_oblique_combinations}));

    end
    i_set = [psyKernel_corr(idx).i_set];
    % for n = 1:numel(x_avg)
    %     text(x_avg(n), y_avg(n), sprintf('%d', i_set(n)),'fontsize',16)
    % end
    plot(x_avg, y_avg, '.','markersize',20, 'color','blue');
    if doErrorbar
        errorbar(x_avg, y_avg,y_sem,y_sem,x_sem,x_sem,'.','color','blue')
    end
    set(gca,'fontsize',16)
    xlabel('Cardinal');
    ylabel('Oblique')
    title([subjectCode, ' corr-model-empirical'])

    
    idx = strcmp({psyKernel_baseline_corr(:).subjectCode}, subjectCode);
    
    switch plotOptions
        case 'corr'
            x_baseline = [psyKernel_baseline_corr(idx).r_model_empirical_cardinal_avg];
            y_baseline = [psyKernel_baseline_corr(idx).r_model_empirical_oblique_avg];
            x_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).r_model_empirical_cardinal_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).r_model_empirical_cardinal_combinations}));

            y_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).r_model_empirical_oblique_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).r_model_empirical_oblique_combinations}));
        case 'dist'
            x_baseline = [psyKernel_baseline_corr(idx).dist_model_empirical_cardinal_avg];
            y_baseline = [psyKernel_baseline_corr(idx).dist_model_empirical_oblique_avg];
            
            x_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).dist_model_empirical_cardinal_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).dist_model_empirical_cardinal_combinations}));

            y_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).dist_model_empirical_oblique_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).dist_model_empirical_oblique_combinations}));
    end
    
    plot(x_baseline, y_baseline, '.','markersize',20, 'color',[0.5, 0.5, 0.5])
    if doErrorbar
        errorbar(x_baseline, y_baseline,y_baseline_sem,y_baseline_sem,x_baseline_sem,x_baseline_sem,'.','color',[0.5, 0.5, 0.5])
    end
    %xlim([0,1]); ylim([0,1])
      dist_baseline   = pdist2([x_baseline;y_baseline]', [0;0]');
    dist_model      = pdist2([x_avg;y_avg]', [0;0]');
    [h,p] = ttest2(dist_baseline,dist_model)

    title(sprintf('%s, %s-model-empirical',plotOptions, subjectCode));

end

%% visualization model-ideal
plotOptions = 'corr'; % corr or dist
doErrobars = true;
figure; 
for i  = 1:2
    subplot(1,2,i); hold on

    subjectCode =  subject_list{i};
    idx = strcmp({psyKernel_corr(:).subjectCode}, subjectCode);
    switch plotOptions
        case 'corr'
            x_avg = [psyKernel_corr(idx).r_model_ideal_cardinal_avg];
            y_avg = [psyKernel_corr(idx).r_model_ideal_oblique_avg];

            x_sem = cellfun(@std,{psyKernel_corr(idx).r_model_ideal_cardinal_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).r_model_ideal_cardinal_combinations}));
            y_sem = cellfun(@std,{psyKernel_corr(idx).r_model_ideal_oblique_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).r_model_ideal_oblique_combinations}));

        case 'dist'
            x_avg = [psyKernel_corr(idx).dist_model_ideal_cardinal_avg];
            y_avg = [psyKernel_corr(idx).dist_model_ideal_oblique_avg];

            x_sem = cellfun(@std,{psyKernel_corr(idx).dist_model_ideal_cardinal_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).dist_model_ideal_cardinal_combinations}));
            y_sem = cellfun(@std,{psyKernel_corr(idx).dist_model_ideal_oblique_combinations}) ./ ...
                    sqrt(cellfun(@numel, {psyKernel_corr(idx).dist_model_ideal_oblique_combinations}));

    end
    % i_set = [psyKernel_corr(idx).i_set];
    % scatter(x_avg, y_avg)
    % for n = 1:numel(x_avg)
    %     text(x_avg(n), y_avg(n), sprintf('%d', i_set(n)),'fontsize',16)
    % end
    plot(x_avg, y_avg, '.','markersize',20, 'color','blue');
    if doErrorbar
        errorbar(x_avg, y_avg,y_sem,y_sem,x_sem,x_sem,'.','color','blue')
    end

    set(gca,'fontsize',16)
    xlabel('Cardinal');
    ylabel('Oblique')
    title([subjectCode, ' corr-model-ideal'])
    
    switch plotOptions
        case 'corr'
            x_empirical_ideal_cardinal_avg = psyKernel_corr(find(idx,1)).r_empirical_ideal_cardinal_avg;
            y_empirical_ideal_oblique_avg  = psyKernel_corr(find(idx,1)).r_empirical_ideal_oblique_avg;


        case 'dist'
             x_empirical_ideal_cardinal_avg = psyKernel_corr(find(idx,1)).dist_empirical_ideal_cardinal_avg;
            y_empirical_ideal_oblique_avg  = psyKernel_corr(find(idx,1)).dist_empirical_ideal_oblique_avg;
    end

            
    plot(x_empirical_ideal_cardinal_avg, y_empirical_ideal_oblique_avg, '*', 'MarkerSize',20,'Color','red');
    if doErrorbar
        errorbar(x_baseline, y_baseline,y_baseline_sem,y_baseline_sem,x_baseline_sem,x_baseline_sem,'.','color',[0.5, 0.5, 0.5])
    end
   % addCircularGrid([r_empirical_ideal_cardinal_avg,r_empirical_ideal_oblique_avg], [0:0.05:1]);
    
    idx = strcmp({psyKernel_baseline_corr(:).subjectCode}, subjectCode);
    


    switch plotOptions
        case 'corr'
            x_baseline = [psyKernel_baseline_corr(idx).r_model_ideal_cardinal_avg];
            y_baseline = [psyKernel_baseline_corr(idx).r_model_ideal_oblique_avg];

            x_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).r_model_ideal_cardinal_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).r_model_ideal_cardinal_combinations}));

            y_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).r_model_ideal_oblique_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).r_model_ideal_oblique_combinations}));
        case 'dist'
            x_baseline = [psyKernel_baseline_corr(idx).dist_model_ideal_cardinal_avg];
            y_baseline = [psyKernel_baseline_corr(idx).dist_model_ideal_oblique_avg];

            x_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).dist_model_ideal_cardinal_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).dist_model_ideal_cardinal_combinations}));

            y_baseline_sem = cellfun(@std,{psyKernel_baseline_corr(idx).dist_model_ideal_oblique_combinations}) ./ ...
                        sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).dist_model_ideal_oblique_combinations}));
    end
    

    
    plot(x_baseline, y_baseline, '.','markersize',20, 'color',[0.5, 0.5, 0.5])
    %xlim([0,1]); ylim([0,1])
    
    %  [h,p] = ttest2(x_baseline, x_avg)
    % [h,p] = ttest2(y_baseline, y_avg)
    dist_baseline   = pdist2([x_baseline;y_baseline]', [r_empirical_ideal_cardinal_avg;r_empirical_ideal_oblique_avg]');
    dist_model      = pdist2([x_avg;y_avg]', [r_empirical_ideal_cardinal_avg;r_empirical_ideal_oblique_avg]');
    [h,p] = ttest2(dist_baseline,dist_model)

    title(sprintf('%s, %s-ideal',plotOptions, subjectCode));
end
%% directly measure assmetry as the distance between kernels

[asymmetry_empirical_avg_Ro, asymmetry_empirical_combinations_Ro] = ...
    compute_distance_kernels(w_ori_cardinal_Ro, w_ori_oblique_Ro, w_ori_cardinal_ideal);
[asymmetry_empirical_avg_Gr, asymmetry_empirical_combinations_Gr] = ...
    compute_distance_kernels(w_ori_cardinal_Gr, w_ori_oblique_Gr, w_ori_cardinal_ideal);
figure
for i  = 1:2
    subplot(1,2,i); hold on

    subjectCode =  subject_list{i};
    idx = strcmp({psyKernel_corr(:).subjectCode}, subjectCode);
    
    assymetry_avg_bestParams = [psyKernel_corr(idx).asymmetry_model_avg];
    assymetry_avg_bestParams_sem = cellfun(@std, {psyKernel_corr(idx).asymmetry_model_combinations}) ./ ...     
                                    sqrt(cellfun(@numel, {psyKernel_corr(idx).asymmetry_model_combinations}));




    i_set = [psyKernel_corr(idx).i_set];
  

    idx = strcmp({psyKernel_baseline_corr(:).subjectCode}, subjectCode);
    assymetry_avg_baseline = [psyKernel_baseline_corr(idx).asymmetry_model_avg];
    assymetry_avg_baseline_sem = cellfun(@std, {psyKernel_baseline_corr(idx).asymmetry_model_combinations}) ./ ...     
                                    sqrt(cellfun(@numel, {psyKernel_baseline_corr(idx).asymmetry_model_combinations}));



    switch subjectCode
        case 'monkeyR'
            asymmetry_empirical_avg     =  asymmetry_empirical_avg_Ro;
        case 'monkeyG'
            asymmetry_empirical_avg     = asymmetry_empirical_avg_Gr;
    end

    assymetry_avg_all   = [assymetry_avg_bestParams, assymetry_avg_baseline];
    is_baseline         = boolean([zeros(size(assymetry_avg_bestParams)), ones(size(assymetry_avg_baseline))]);
    assymetry_sem_all   = [assymetry_avg_bestParams_sem, assymetry_avg_baseline_sem];
    
    [~, idx_sort] = sort(assymetry_avg_all, 'ascend');
    
    assymetry_avg_sorted    = assymetry_avg_all(idx_sort);
    assymetry_sem_sorted    = assymetry_sem_all(idx_sort);
    is_baseline_sorted      = is_baseline(idx_sort);

    plot(find(is_baseline_sorted), assymetry_avg_sorted(is_baseline_sorted),'.','markersize',20, 'color',[0.5, 0.5, 0.5]); hold on
    plot(find(~is_baseline_sorted), assymetry_avg_sorted(~is_baseline_sorted),'.','markersize',20, 'color','blue');
    yl  =yline(asymmetry_empirical_avg,'--','empirical value');
    yl.FontSize = 1; 
    if doErrorbar
        errorbar(find(is_baseline_sorted), assymetry_avg_sorted(is_baseline_sorted), assymetry_sem_sorted(is_baseline_sorted),'.', ...
            'color',[0.5, 0.5, 0.5]);
        errorbar(find(~is_baseline_sorted), assymetry_avg_sorted(~is_baseline_sorted), assymetry_sem_sorted(~is_baseline_sorted),'.', ...
            'color','blue');
    end
    
    set(gca,'fontsize',18);
    ylabel('Cardinal Kernel - Oblique Kernel')
    title(sprintf('%s, kernel asymmetry',subjectCode));


end

%%
figure
%for i = 1:2
% subplot(1,2,i); hold on
%
% subjectCode =  subject_list{i};
% idx = strcmp({psyKernel_corr(:).subjectCode}, subjectCode);


asymmetry_prior_bestParams  = [psyKernel_corr(:).asymmetry_prior];
asymmetry_delta_bestParams  = [psyKernel_corr(:).asymmetry_delta];
prior_cardinal_bestParams   = [psyKernel_corr(:).prior_cardinal];
prior_oblique_bestParams    = [psyKernel_corr(:).prior_oblique];
delta_cardinal_bestParams   = [psyKernel_corr(:).delta_cardinal];
delta_oblique_bestParams    = [psyKernel_corr(:).delta_oblique];
assymetry_avg_bestParams    = [psyKernel_corr(:).asymmetry_model_avg];

min_prior_bestParams        = min([prior_cardinal_bestParams;prior_oblique_bestParams],[],1);
min_delta_bestParams        =  min([delta_cardinal_bestParams ; delta_oblique_bestParams],[],1);

% i_set = [psyKernel_corr(idx).i_set];


%idx = strcmp({psyKernel_baseline_corr(:).subjectCode}, subjectCode);
assymetry_avg_baseline = [psyKernel_baseline_corr(:).asymmetry_model_avg];
asymmetry_prior_baseline  = [psyKernel_baseline_corr(:).asymmetry_prior];
asymmetry_delta_baseline  = [psyKernel_baseline_corr(:).asymmetry_delta];
prior_cardinal_baseline   = [psyKernel_baseline_corr(:).prior_cardinal];
prior_oblique_baseline    = [psyKernel_baseline_corr(:).prior_oblique];
delta_cardinal_baseline   = [psyKernel_baseline_corr(:).delta_cardinal];
delta_oblique_baseline    = [psyKernel_baseline_corr(:).delta_oblique];

min_prior_baseline        = min([prior_cardinal_baseline;prior_oblique_baseline],[],1);
min_delta_baseline        =  min([delta_cardinal_baseline ; delta_oblique_baseline],[],1);

% switch subjectCode
%     case 'monkeyR'
%         asymmetry_empirical_avg     =  asymmetry_empirical_avg_Ro;
%     case 'monkeyG'
%         asymmetry_empirical_avg     = asymmetry_empirical_avg_Gr;
% end

%%%% two dimensional heatmap

%
x_list = {'asymmetry_prior';'asymmetry_delta';'min_prior';'min_delta'};%;'delta_cardinal';'delta_oblique'};
figure
for i = 1:4
subplot(2,2,i); hold on
eval(sprintf('x_bestParams = %s_bestParams;',x_list{i}));
eval(sprintf('x_baseline = %s_baseline;',x_list{i}));

y_bestParams = assymetry_avg_bestParams;
y_baseline = assymetry_avg_baseline;
if ismember(x_list{i}, {'min_prior';'min_delta'})
    y_bestParams = abs(y_bestParams);
    y_baseline = abs(y_baseline);
end
plot(x_bestParams, y_bestParams,'.','markersize',20,'Color','blue');

plot(x_baseline, y_baseline,'.','markersize',20,'Color',[0.5, 0.5, 0.5]);
xlabel(x_list{i},'interpreter','none'); ylabel('asymmetry_kernel','interpreter','none')
set(gca,'fontsize',18);

x_all = [x_bestParams, x_baseline];
y_all = [y_bestParams, y_baseline];


p = polyfit(x_all, y_all, 1);

% 3. Evaluate the polynomial over a smoother range of points for the curve
x_fit = linspace(min(x_all), max(x_all), 200);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '--', 'LineWidth', 2);
end






%end
%% visualization of the model kernels
doThis = 1;
if doThis
    figure;
    %%% best params monkey R
     subjectCode = 'monkeyR';
    idx = find(strcmp({psyKernel_corr(:).subjectCode}, subjectCode));
    for n = 1:10
        subplot(3,4,n); hold on
        plot_kernel_simple(psyKernel_corr, idx(n));
    end
    sgtitle('Best params, monkey R')


    figure
    %%%% best params monkey G
    subjectCode = 'monkeyG';
   
    idx = find(strcmp({psyKernel_corr(:).subjectCode}, subjectCode));
    for n = 1:10
        subplot(3,4,n); hold on
        plot_kernel_simple(psyKernel_corr, idx(n));
    end
    sgtitle('Best params, monkey G')





    figure;
    subjectCode = 'monkeyG';
   
    idx = find(strcmp({psyKernel_baseline_corr(:).subjectCode}, subjectCode));
    for n = 1:30
        subplot(5,6,n); hold on
        plot_kernel_simple(psyKernel_baseline_corr, idx(n));
       
    end
    sgtitle('Symmetric params')
end
%% symmetry as a function of prior, delta

%% check how good these selected sets are in terms of matching neural data
doThis = 0;
if doThis

    load('../../results/model_match_neural_rolo_moreRepeats.mat');
    load('../../results/model_match_neural_gremlin_moreRepeats.mat');
    
    [~, index_rolo] = sort(neural_match_result_rolo.dist_feature, 'ascend');
    [~, index_gremlin] = sort(neural_match_result_gremlin.dist_feature, 'ascend');
    
    monkeyR_set_list = [7,8];
    monkeyG_set_list = [4,8];
    nPlot = 30;
    
    
    
    %%% symmetric explain neural?
    nSet = numel(psyKernel_baseline_corr)/ 2;
    [idx_baseline_rolo, idx_baseline_gremlin, dist_baseline_rolo, dist_baseline_gremlin] = deal(zeros(nSet,1));
    for i = 1:nSet
        prior_cardinal = psyKernel_baseline_corr(i).prior_cardinal;
        prior_oblique = psyKernel_baseline_corr(i).prior_oblique;
        delta_cardinal = psyKernel_baseline_corr(i).delta_cardinal;
        delta_oblique = psyKernel_baseline_corr(i).delta_oblique;
    
        idx = find(neural_match_result_rolo.prior_cardinal == prior_cardinal & ...
              neural_match_result_rolo.prior_oblique == prior_oblique & ...
              neural_match_result_rolo.delta_cardinal == delta_cardinal & ...
              neural_match_result_rolo.delta_oblique  == delta_oblique);
    
        dist_baseline_rolo(i) = neural_match_result_rolo.dist_feature(idx);
        idx_baseline_rolo(i) = find(index_rolo == idx);
    
        idx = find(neural_match_result_gremlin.prior_cardinal == prior_cardinal & ...
              neural_match_result_gremlin.prior_oblique == prior_oblique & ...
              neural_match_result_gremlin.delta_cardinal == delta_cardinal & ...
              neural_match_result_gremlin.delta_oblique  == delta_oblique);
    
        dist_baseline_gremlin(i) = neural_match_result_gremlin.dist_feature(idx);
        idx_baseline_gremlin(i) = find(index_gremlin == idx);
    
    end
    
    
    figure;
    subplot(2,3,[1,2]); hold on
    plot(neural_match_result_rolo.dist_feature(index_rolo),'LineWidth',2);
    plot(idx_baseline_rolo, dist_baseline_rolo, '.','markersize',20, 'color',[0.5, 0.5, 0.5])
    set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); xlim([0,900])
    title('Monkey R')
    subplot(2,3,3); hold on
    plot(neural_match_result_rolo.dist_feature(index_rolo(1:nPlot)),'-o');
    for n = 1:numel(monkeyR_set_list)
        i = monkeyR_set_list(n);
        plot(i, neural_match_result_rolo.dist_feature(index_rolo(i)),'*','color','red');
    end
    set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); 
    
    
    subplot(2,3,[4,5]); hold on
    plot(neural_match_result_gremlin.dist_feature(index_gremlin),'LineWidth',2);
    plot(idx_baseline_gremlin, dist_baseline_gremlin, '.','markersize',20, 'color',[0.5, 0.5, 0.5])
    set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); xlim([0,900])
    title('Monkey G')
    subplot(2,3,6); hold on
    plot(neural_match_result_gremlin.dist_feature(index_gremlin(1:nPlot)),'-o');
    for n = 1:numel(monkeyR_set_list)
        i = monkeyR_set_list(n);
        plot(i, neural_match_result_gremlin.dist_feature(index_gremlin(i)),'*','color','red');
    end
    set(gca,'fontsize',18); xlabel('Index of param set'); ylabel('Distance'); 
end
%% helper functions
function [r_avgs, r_combinations] = compute_corr_kernels(w_ori_1, w_ori_2)
if size(w_ori_1,1) ~= size(w_ori_2,1)
    error('Size of psychometric kernels need to match')
end
w_ori_1_avg = mean(w_ori_1, 2);
w_ori_2_avg = mean(w_ori_2, 2);

r_avgs = corr(w_ori_1_avg, w_ori_2_avg);

N_1 = size(w_ori_1, 2);
N_2 = size(w_ori_2, 2);
r_combinations = zeros(N_1, N_2);
for n1 = 1:N_1
    for n2 = 1:N_2
        r_combinations(n1, n2) = corr(w_ori_1(:, n1), w_ori_2(:, n2)); 
    end
end
r_combinations = r_combinations(:);

end

function [dist_avgs, dist_combinations] = compute_distance_kernels(w_ori_cardinal, w_ori_oblique, w_ori_cardinal_ideal)

if size(w_ori_cardinal,1) ~= size(w_ori_oblique,1)
    error('Size of psychometric kernels need to match')
end
%%% shift the oblique kernel by 45 degrees
w_ori_oblique_shifted =  circshift(w_ori_oblique, -3, 1);
%%% flip sign based on cardinal_ideal
w_ori_cardinal_flipped = w_ori_cardinal .* sign(w_ori_cardinal_ideal);
w_ori_oblique_flipped = w_ori_oblique_shifted .* sign(w_ori_cardinal_ideal);

% 
w_ori_cardinal_avg = mean(w_ori_cardinal_flipped, 2);
w_ori_oblique_avg = mean(w_ori_oblique_flipped, 2);

%w_cardinal_avg =

dist_avgs = sum(w_ori_cardinal_avg - w_ori_oblique_avg);



N_1 = size(w_ori_cardinal_flipped, 2);
N_2 = size(w_ori_oblique_flipped, 2);
dist_combinations = zeros(N_1, N_2);
for n1 = 1:N_1
    for n2 = 1:N_2
        dist_combinations(n1, n2) = sum(w_ori_cardinal_flipped(:, n1) - w_ori_oblique_flipped(:, n2));
    end
end
dist_combinations = dist_combinations(:);

%%% normalize by an ideal kernel and zero
dist_norm = sum(abs(w_ori_cardinal_ideal));

dist_avgs = dist_avgs / dist_norm;
dist_combinations = dist_combinations / dist_norm;
end



function plot_psykernel(subjectCode, i_set, kernel_type)
data_folder = '../../results/model_behavior/psyKernel_bestParams';
data_name = fullfile(data_folder, sprintf('psyKernels_best_parameters_%s_set_%d', subjectCode, i_set));
load(data_name);
for n = 1:numel(psyKernel_model)
    %%%% add these fields just for compaibility of the plot function
    psyKernel_model(n).mode = 'Single';
    psyKernel_model(n).image_task = psyKernel_model(n).task;
end
idx_cardinal    = strcmp({psyKernel_model(:).task},'cardinal');
idx_oblique     = strcmp({psyKernel_model(:).task},'oblique');

prior_cardinal_plot = psyKernel_model(find(idx_cardinal,1)).prior;
delta_cardinal_plot = psyKernel_model(find(idx_cardinal,1)).delta;

prior_oblique_plot = psyKernel_model(find(idx_oblique,1)).prior;
delta_oblique_plot = psyKernel_model(find(idx_oblique,1)).delta;


% figure
% subplot(1,2,1)
plotOptions.mode            = 'Single';
plotOptions.image_task      = 'cardinal';
plotOptions.prior_task      = prior_cardinal_plot;
plotOptions.delta           = delta_cardinal_plot;
plotOptions.plotIndividual  = false;
fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, plotOptions);
% title(sprintf('Delta = %.2f,  Prior = %.1f', ...
%     plotOptions.delta, plotOptions.prior_task));

%subplot(1,2,2)
plotOptions.mode            = 'Single';
plotOptions.image_task      = 'oblique';
plotOptions.prior_task      = prior_oblique_plot;
plotOptions.delta           = delta_oblique_plot;
plotOptions.plotIndividual  = false;
fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, plotOptions);
% title(sprintf('Delta = %.2f,  Prior = %.1f', ...
%     plotOptions.delta, plotOptions.prior_task));

title({sprintf('prior_{cardinal} = %.1f, delta_{cardinal} = %.2f', prior_cardinal_plot, delta_cardinal_plot);...
       sprintf('prior_{oblique} = %.1f, delta_{oblique} = %.2f', prior_oblique_plot, delta_oblique_plot) })
% sgtitle(sprintf('Parameters set %d for %s',i_set, subjectCode),'fontsize',20,'fontweight','bold')
end



function addCircularGrid(center, radii)

x0 = center(1);
y0 = center(2);

theta = linspace(0,2*pi,300);

hold on

for r = radii
    x = x0 + r*cos(theta);
    y = y0 + r*sin(theta);
    plot(x,y,'Color','k','LineStyle','--','LineWidth',1)
end

end

function plot_kernel_simple(psyKernel_corr, i)

w_cardinal_model =  psyKernel_corr(i).w_cardinal_model;
w_oblique_model =   psyKernel_corr(i).w_oblique_model;

w_cardinal_empirical =    psyKernel_corr(i).w_cardinal_empirical;
w_oblique_empirical =   psyKernel_corr(i).w_oblique_empirical;

w_cardinal_ideal =    psyKernel_corr(i).w_cardinal_ideal;
w_oblique_ideal=     psyKernel_corr(i).w_oblique_ideal;


plot(w_cardinal_model,'color','red','LineWidth',2);
plot(w_oblique_model,'color','blue','LineWidth',2);

plot(w_cardinal_empirical,'color','red','linestyle','--','LineWidth',2);
plot(w_oblique_empirical,'color','blue','linestyle','--','LineWidth',2);

plot(w_cardinal_ideal,'color','black','linestyle','--','LineWidth',2);
plot(w_oblique_ideal,'color','black','linestyle','--','LineWidth',2);
title(sprintf('asym = %.2f',psyKernel_corr(i).asymmetry_model_avg))
ylim([-0.5, 0.5])
end