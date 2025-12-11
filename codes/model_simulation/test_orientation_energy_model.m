clear all
clc
close all
%%
nX = 32;
P = S_Exp_Para('test-interleaved','G.fct','nxN','G.nx',nX);
projective_fields = C_Projection(P.G.fct, P.G.nx, P.G.dimension_X, P.G.dimension_G, P.G.number_locations);

im_type     = P.I.fct;
n_locs      = P.G.number_locations;
im_height   = P.G.ny;


kernelOptions.orientationsDEG      = [0:15:165]; 
kernelOptions.filterMode           = 'hard';
kernelOptions.filterSDDeg          = 15;


%% 1. four orientations, difference contrast
nTrial                  = 50;
num_ori_bin             = 12;

figure
subplot(2,2,1)
image_task               = 'cardinal';
stimulus_contrast_list  = [0:5:15];
plot_ori_energy(image_task, stimulus_contrast_list,nTrial, kernelOptions, im_type, n_locs, im_height)

subplot(2,2,2)
image_task               = 'cardinal';
stimulus_contrast_list   = -[0:5:15];
plot_ori_energy(image_task, stimulus_contrast_list,nTrial, kernelOptions, im_type, n_locs, im_height)

subplot(2,2,3)
image_task               = 'oblique';
stimulus_contrast_list  = [0:5:15];
plot_ori_energy(image_task, stimulus_contrast_list,nTrial, kernelOptions, im_type, n_locs, im_height)

subplot(2,2,4)
image_task               = 'oblique';
stimulus_contrast_list   = -[0:5:15];
plot_ori_energy(image_task, stimulus_contrast_list,nTrial, kernelOptions, im_type, n_locs, im_height)




%%
function plot_ori_energy(image_task, stimulus_contrast_list,nTrial, kernelOptions, im_type, n_locs, im_height)
num_ori_bin = numel(kernelOptions.orientationsDEG);
for n = 1:numel(stimulus_contrast_list)
    orientation_energy = zeros(nTrial, num_ori_bin);
    if stimulus_contrast_list(n) < 0
        stimulus_contrast = [stimulus_contrast_list(n), 0];
    elseif stimulus_contrast_list(n) == 0
        stimulus_contrast = [0,0];
    else
        stimulus_contrast = [0, stimulus_contrast_list(n)];
    end

    for t = 1:nTrial
        stim_image = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);

        % stimuliF = fftshift(fft2(stim_image));
        % 
        % for iOri = 1:num_ori_bin
        %     orientation_energy(t,iOri) = ...
        %         computesignal(stimuliF, kernelOptions.filterMode, ...
        %         kernelOptions.orientationsDEG(iOri), kernelOptions.filterSDDeg);
        % end

        [orientation_bins_center, orientation_energy(t,:)] = get_ori_energy_model(stim_image, num_ori_bin);

    end
    %orientation_bins = kernelOptions.orientationsDEG;
    orientation_bins = orientation_bins_center / pi * 180;
    avg = mean(orientation_energy, 1);
    sdv = std(orientation_energy, [], 1);
    h = plot(orientation_bins,avg, 'LineWidth', 2); hold on

    fill([orientation_bins, fliplr(orientation_bins)],...
        [avg - sdv, fliplr(avg + sdv)],h.Color,'FaceAlpha',0.5)

end

end