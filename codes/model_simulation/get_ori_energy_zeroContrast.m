clear all
clc
close all
%%%% Generate zero-contrast images for normalization
%%
nTrials = 1000;
nX      = 32;
P       = S_Exp_Para('test-interleaved','G.fct','nxN');

nNeuron     = P.G.dimension_X; 
im_type     = P.I.fct;
n_locs      = P.G.number_locations;
im_height   = P.G.ny;

stimulus_contrast   = [0,0];
image_task          = 'cardinal';

num_ori_bin         = 12;
phi_x               = P.G.phi_x;
%%
n_pixels = size(P.G.G,1);
zero_energy_proj = zeros(nTrials, nNeuron);
zero_energy_fft  = zeros(nTrials, num_ori_bin);
for t = 1:nTrials
    stim_image = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);

    %%%% proj
    proj = (stim_image(:)' * P.G.G) / n_pixels;
    zero_energy_proj(t,:) = proj .^ 2;
    %%%% fft
    [orientation_bins_center, zero_energy_fft(t,:)] = get_ori_energy_model(stim_image, num_ori_bin);

end

orientation_bins_center = 180 * orientation_bins_center / pi;
phi_x                   = 180 * phi_x / pi; 

save('../../results/model_oriEnergy_zeroContrast',...
    'zero_energy_proj','zero_energy_fft','orientation_bins_center','phi_x');