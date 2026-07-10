
clear all
clc
close all
%% generate projective fields (I only need the phix from this step)
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
phi_x = 180 * phi_x / pi; % convert from [0, pi] to [0, 180] 
%% 
cardinal_1 = 0;
cardinal_2 = 90;
dist_to_cardinal_1 = abs(mod((phi_x - cardinal_1) + 90, 180) - 90);
dist_to_cardinal_2 = abs(mod((phi_x - cardinal_2) + 90, 180) - 90);
dist_to_cardinal = min([dist_to_cardinal_1;dist_to_cardinal_2],[],1);


oblique_1 = 45;
oblique_2 = 135;
dist_to_oblique_1 = abs(mod((phi_x - oblique_1) + 90, 180) - 90);
dist_to_oblique_2 = abs(mod((phi_x - oblique_2) + 90, 180) - 90);
dist_to_oblique = min([dist_to_oblique_1;dist_to_oblique_2],[],1);

%%
highC_highO  = find(dist_to_cardinal < prctile(dist_to_cardinal, 25) & ...
                 dist_to_oblique < prctile(dist_to_oblique, 25));

lowC_lowO  = find(dist_to_cardinal >= prctile(dist_to_cardinal, 75) & ...
                 dist_to_oblique >= prctile(dist_to_oblique, 75));

