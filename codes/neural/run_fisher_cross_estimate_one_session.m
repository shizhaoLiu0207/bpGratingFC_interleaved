function [dat_cross_fisher,small_trial_flag] = run_fisher_cross_estimate_one_session(dat_cross_fisher,data_run,info_run)
%%% This function runs bias-corrected estimation of FI for one session and saves results
%%% in a big struct.
%%% Within this session, run separately for each coherence level and each
%%% task (cardinal/oblique)
%%% Shizhao Liu 10/09/2024


signal              = data_run.signal;
orientation         = data_run.orientation;
spikeCount          = data_run.spikeCount;
sessionStr          = info_run.sessionStr;
sessionType         = info_run.sessionType;
timeWin             = info_run.timeWin;
i_win               = info_run.i_win;


signal_abs_list     = nonzeros(unique(abs(signal)));
nNeuron             = size(spikeCount,2);

for i = 1:numel(signal_abs_list) % run analysis for each coherence level
    coher               = signal_abs_list(i);

    % separate by task
    idx_cardinal_pos    = find(signal == coher & ismember(orientation,[0,90]));
    idx_cardinal_neg    = find(signal == -coher & ismember(orientation,[0,90]));
    idx_oblique_pos     = find(signal == coher & ismember(orientation,[45,135]));
    idx_oblique_neg     = find(signal == -coher & ismember(orientation,[45,135]));

    

    idx_cardinal_zero = find(signal == 0 & ismember(orientation,[0,90]));
    idx_oblique_zero  = find(signal == 0 & ismember(orientation,[45,135]));

    %%% zero coherence activity
    T_zero_cardinal   = numel(idx_cardinal_zero);
    T_zero_oblique    = numel(idx_oblique_zero); 

    X_zero_cardinal   = spikeCount(idx_cardinal_zero,:);
    X_zero_oblique    = spikeCount(idx_oblique_zero,:); 

    %%%%%%% non-coherence activity
    T1_cardinal = numel(idx_cardinal_pos);
    T2_cardinal = numel(idx_cardinal_neg);
    T1_oblique = numel(idx_oblique_pos);
    T2_oblique = numel(idx_oblique_neg);


    X_cardinal_pos_all           = spikeCount(idx_cardinal_pos,:); 
    X_cardinal_neg_all           = spikeCount(idx_cardinal_neg,:);
    X_cardinal_all               = [X_cardinal_pos_all;X_cardinal_neg_all];
    Y_cardinal_all               = [ones(T1_cardinal,1); -1 * ones(T2_cardinal,1)];

    X_oblique_pos_all           = spikeCount(idx_oblique_pos,:); 
    X_oblique_neg_all           = spikeCount(idx_oblique_neg,:);
    X_oblique_all               = [X_oblique_pos_all;X_oblique_neg_all];
    Y_oblique_all               = [ones(T1_oblique,1); -1 * ones(T2_oblique,1)];

    
    %%%%%% cardinal task fisher information, using Cov from cardinal zero cohernece trials %%%%%%%%%%
    
   % results_crossI = bias_corrected_fisherInfo_cross(X_cardinal_all, Y_cardinal_all, X_zero_cardinal);

    if (T_zero_cardinal < (nNeuron + 4)) | (T_zero_oblique < (nNeuron + 4)) 
       
        continue
        
    end

    
    results_crossI = bias_corrected_fisherInfo_cross(X_cardinal_all, Y_cardinal_all, ...
                                                            X_oblique_all, Y_oblique_all,...
                                                            X_zero_cardinal, X_zero_oblique, coher);


    i_d = numel(dat_cross_fisher) + 1;
    dat_cross_fisher(i_d).sessionStr            = sessionStr;
    dat_cross_fisher(i_d).sessionType           = sessionType;
    dat_cross_fisher(i_d).nNeuron               = nNeuron;
    dat_cross_fisher(i_d).coherence_level       = coher;
    dat_cross_fisher(i_d).task                  = 'interleaved';
  

    dat_cross_fisher(i_d).T1_cardinal           = T1_cardinal;
    dat_cross_fisher(i_d).T2_cardinal           = T2_cardinal;

    dat_cross_fisher(i_d).T1_oblique            = T1_oblique;
    dat_cross_fisher(i_d).T2_oblique            = T2_oblique;

    dat_cross_fisher(i_d).T_zero_cardinal       =  T_zero_cardinal;
    dat_cross_fisher(i_d).T_zero_oblique        =  T_zero_oblique;

    dat_cross_fisher(i_d).N                     = nNeuron;

    dat_cross_fisher(i_d).timeWin               = timeWin;
    dat_cross_fisher(i_d).timeWinIndex          = i_win - 1;  

    dat_cross_fisher(i_d).fisher_cardinal_cardinal_bc         = results_crossI.I_AA_bc;
    dat_cross_fisher(i_d).fisher_oblique_oblique_bc         = results_crossI.I_BB_bc;

    dat_cross_fisher(i_d).fisher_cardinal_oblique_bc         = results_crossI.I_AB_bc;
    dat_cross_fisher(i_d).fisher_oblique_cardinal_bc         = results_crossI.I_BA_bc;
   

    dat_cross_fisher(i_d).fisher_cardinal_cardinal_shuffle_bc       = results_crossI.I_AA_shuffle_bc;
    dat_cross_fisher(i_d).fisher_oblique_oblique_shuffle_bc         = results_crossI.I_BB_shuffle_bc;

    dat_cross_fisher(i_d).fisher_cardinal_oblique_shuffle_bc         = results_crossI.I_AB_shuffle_bc;
    dat_cross_fisher(i_d).fisher_oblique_cardinal_shuffle_bc         = results_crossI.I_BA_shuffle_bc;
   


    dat_cross_fisher(i_d).var_cardinal_cardinal_bc           = results_crossI.var_IAA_bc;
    dat_cross_fisher(i_d).var_oblique_oblique_bc             = results_crossI.var_IBB_bc;

    dat_cross_fisher(i_d).var_cardinal_oblique_bc             = results_crossI.var_IAB_bc;
    dat_cross_fisher(i_d).var_oblique_cardinal_bc             = results_crossI.var_IBA_bc;


    dat_cross_fisher(i_d).var_cardinal_cardinal_shuffle_bc             = results_crossI.var_IAA_bc_shuffle;
    dat_cross_fisher(i_d).var_oblique_oblique_shuffle_bc             = results_crossI.var_IBB_bc_shuffle;

    dat_cross_fisher(i_d).var_cardinal_oblique_shuffle_bc             = results_crossI.var_IAB_bc_shuffle;
    dat_cross_fisher(i_d).var_oblique_cardinal_shuffle_bc             = results_crossI.var_IBA_bc_shuffle;

   
   
end
end


 