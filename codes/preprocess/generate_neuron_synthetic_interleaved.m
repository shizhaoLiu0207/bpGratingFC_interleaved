
clear all
clc
close all
%% generate projective fields

%%%%% Only b_PF and nNeuron are useful parameters
b_PF_original        = 0; % 0 or 0.8
nNeuron     = 256;
%%%%%%% these task parameters are not used since we don't need the model to
%%%%%%% perform the task. Specifying them just for no errors
image_task              = 'cardinal';
prior_task              = [0,1];
delta                   = 0.08;
stimulus_contrast       = [0,0];
              
P = S_Exp_Para('run-interleaved-nonuniform', 'I.stimulus_contrast',stimulus_contrast,...
    'G.prior_task',prior_task,...
    'I.image_task',image_task,...
    'G.b_PF',b_PF_original,...
    'G.delta',delta,...
    'G.dimension_X',nNeuron);

projective_fields = C_Projection(P.G.fct, P.G.nx, P.G.dimension_X, P.G.dimension_G, P.G.number_locations, P.G.b_PF);


phi_x = projective_fields.phi_x; %%% preferred orientation of synthetic neurons

%%
%%%% Scenerio 1: b_PF = 2, tao_cardinal = 0, tao_oblique = 0
%%%% Scenerio 2: b_PF = -2, tao_cardinal = 0, tao_oblique = 0
%%%% Scenerio 3: b_PF = 0, tao_cardinal = 2, tao_oblique = 0
%%%% Scenerio 4: b_PF = 0, tao_cardinal = -2, tao_oblique = 0
%%%% Scenerio 5: b_PF = 0, tao_cardinal = 2, tao_oblique = 2
%%%% Scenerio 6: b_PF = 0, tao_cardinal = 2, tao_oblique = -2
nSample             = 128;

params_list(1).b_PF = 2;    params_list(1).tao_cardinal = 0;    params_list(1).tao_oblique = 0;
params_list(2).b_PF = -2;   params_list(2).tao_cardinal = 0;    params_list(2).tao_oblique = 0;
params_list(3).b_PF = 0;    params_list(3).tao_cardinal = 2;    params_list(3).tao_oblique = 0;
params_list(4).b_PF = 0;    params_list(4).tao_cardinal = -2;   params_list(4).tao_oblique = 0;
params_list(5).b_PF = 0;    params_list(5).tao_cardinal = 2;    params_list(5).tao_oblique = 2;
params_list(6).b_PF = 0;    params_list(6).tao_cardinal = 2;    params_list(6).tao_oblique = -2;

save_folder         = '../../results/filter_neuron_synthetic/subsample';
for k = 1:numel(params_list)
    b_PF_sample     = params_list(k).b_PF;
    tao_cardinal    =  params_list(k). tao_cardinal;
    tao_oblique     = params_list(k).tao_oblique;

    [idx_sample, y, y_c,y_o, y_all] = util_it.sampling_vonMises(phi_x, b_PF_sample, tao_cardinal,...
                                            tao_oblique, nSample);
    phi_x_sample = phi_x(idx_sample);


    %[idx_sample, y, y_c, y_o, y_all] = util_it.sampling_vonMises(theta, b_PF, tao_cardinal, tao_oblique, nSample);
    x = phi_x / pi * 180;
    
    %figure;
   % subplot(2,1,1); hold on 
    % plot(x, y); plot(x, y_c); plot(x, y_o); plot(x, y_all);
    % legend('y','y_c','y_o','y_{all}')
    % xlabel('Orientation (degree)'); ylabel('Probability')
    % set(gca,'fontsize',18);
    
    
    subplot(2,3,k); hold on
    histogram(x,[0:10:180]);
    histogram(x(idx_sample) , [0:10:180]);
    legend('Original','Sampled')
    xlabel('Orientation (degree)');  ylabel('Number of neurons')
    set(gca,'fontsize',18);
    
    title(sprintf('b_{PF,sample} = %.1f, \\tau_{cardinal} = %.1f, \\tau_{oblique} = %.1f',...
        b_PF_sample, tao_cardinal, tao_oblique),'fontsize',20);

    save_name = sprintf('nNeuron_%d_bPF_original_%d_bPF_sample_%d_tao_cardinal_%d_tao_oblique_%d_nSample_%d', ...
                nNeuron, b_PF_original, params_list(k).b_PF, params_list(k).tao_cardinal, params_list(k).tao_oblique, nSample);

    save(fullfile(save_folder,save_name),'idx_sample','phi_x_sample');
end