function predResults = run_cross_predict_model(dat_cardinal,dat_oblique, nFold)

energy_use = 'proj';
kernelOptions.orientationsDEG = [0:15:165]; % orientation bin: 0 to 180, 15 degree increment
kernelOptions.filterMode ='hard';
kernelOptions.filterSDDeg = 15;
fittingKernelParams.hypers  = [0.5 6 60 0.15 0.06];
fittingKernelParams.nLapse  = 2;
[X_cardinal_norm, choice_cardinal, contrast_cardinal] = util_it.extract_X_psyKernel(dat_cardinal.data_use, energy_use);
[X_oblique_norm,  choice_oblique, contrast_oblique]  = util_it.extract_X_psyKernel(dat_oblique.data_use, energy_use);

X.cardinal            = X_cardinal_norm;
X.oblique             = X_oblique_norm;
choice.cardinal       = choice_cardinal;
choice.oblique        = choice_oblique;
choice_list.cardinal  = unique(choice_cardinal);
choice_list.oblique   = unique(choice_oblique);
signal_level.cardinal = contrast_cardinal;
signal_level.oblique  = contrast_oblique;

predResults = psyKernel.fitSPTP_kernel_crosstest(X,choice,kernelOptions,choice_list,signal_level,nFold,fittingKernelParams);

end