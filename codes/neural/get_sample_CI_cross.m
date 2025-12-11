function results_sample = get_sample_CI_cross(results_sample)
nData = numel(results_sample);
CI_level = 68;
prc = [(100-CI_level)/2, 100 - (100-CI_level)/2];

if isfield(results_sample(1), 'combine_fisher_cardinal_cardinal_sample')
    %%%%%% deal with data of combined coherence %%%% 
    for n = 1:nData
        %%%%% median and CI of 4 * 2 fisher information

        results_sample(n).fisher_cardinal_cardinal_median = median(results_sample(n).combine_fisher_cardinal_cardinal_sample, 'omitnan');
        results_sample(n).fisher_oblique_cardinal_median = median(results_sample(n).combine_fisher_oblique_cardinal_sample, 'omitnan');
        results_sample(n).fisher_oblique_oblique_median = median(results_sample(n).combine_fisher_oblique_oblique_sample, 'omitnan');
        results_sample(n).fisher_cardinal_oblique_median = median(results_sample(n).combine_fisher_cardinal_oblique_sample, 'omitnan');
    
        results_sample(n).fisher_cardinal_cardinal_shuffle_median = median(results_sample(n).combine_fisher_cardinal_cardinal_shuffle_sample, 'omitnan');
        results_sample(n).fisher_oblique_cardinal_shuffle_median = median(results_sample(n).combine_fisher_oblique_cardinal_shuffle_sample, 'omitnan');
        results_sample(n).fisher_oblique_oblique_shuffle_median = median(results_sample(n).combine_fisher_oblique_oblique_shuffle_sample, 'omitnan');
        results_sample(n).fisher_cardinal_oblique_shuffle_median = median(results_sample(n).combine_fisher_cardinal_oblique_shuffle_sample, 'omitnan');


        results_sample(n).fisher_cardinal_cardinal_CI = prctile(results_sample(n).combine_fisher_cardinal_cardinal_sample,  prc);
        results_sample(n).fisher_oblique_cardinal_CI = prctile(results_sample(n).combine_fisher_oblique_cardinal_sample,  prc);
        results_sample(n).fisher_oblique_oblique_CI = prctile(results_sample(n).combine_fisher_oblique_oblique_sample,  prc);
        results_sample(n).fisher_cardinal_oblique_CI = prctile(results_sample(n).combine_fisher_cardinal_oblique_sample,  prc);
    
        results_sample(n).fisher_cardinal_cardinal_shuffle_CI = prctile(results_sample(n).combine_fisher_cardinal_cardinal_shuffle_sample,  prc);
        results_sample(n).fisher_oblique_cardinal_shuffle_CI = prctile(results_sample(n).combine_fisher_oblique_cardinal_shuffle_sample,  prc);
        results_sample(n).fisher_oblique_oblique_shuffle_CI = prctile(results_sample(n).combine_fisher_oblique_oblique_shuffle_sample,  prc);
        results_sample(n).fisher_cardinal_oblique_shuffle_CI = prctile(results_sample(n).combine_fisher_cardinal_oblique_shuffle_sample,  prc);
    
    
        %%%% diff between I_cross and I_real as a percentage of I_cross
        
        results_sample(n).combine_diff_info_percent_cardinal_sample = 100 * (results_sample(n).combine_fisher_cardinal_oblique_sample - ...
                results_sample(n).combine_fisher_cardinal_cardinal_sample) ./ results_sample(n).combine_fisher_cardinal_oblique_sample;
        results_sample(n).combine_diff_info_percent_oblique_sample = 100 * (results_sample(n).combine_fisher_oblique_cardinal_sample - ...
                results_sample(n).combine_fisher_oblique_oblique_sample) ./ results_sample(n).combine_fisher_oblique_cardinal_sample;

        results_sample(n).diff_info_percent_cardinal_median = median(results_sample(n).combine_diff_info_percent_cardinal_sample);
        results_sample(n).diff_info_percent_oblique_median  = median(results_sample(n).combine_diff_info_percent_oblique_sample);

        results_sample(n).diff_info_percent_cardinal_CI     = prctile(results_sample(n).combine_diff_info_percent_cardinal_sample, prc); 
        results_sample(n).diff_info_percent_oblique_CI      = prctile(results_sample(n).combine_diff_info_percent_oblique_sample, prc); 

        %%% also do this for the shuffled information
        results_sample(n).combine_diff_shuffle_percent_cardinal_sample = 100 * (results_sample(n).combine_fisher_cardinal_oblique_shuffle_sample - ...
                results_sample(n).combine_fisher_cardinal_cardinal_shuffle_sample) ./ results_sample(n).combine_fisher_cardinal_oblique_shuffle_sample;
        results_sample(n).combine_diff_shuffle_percent_oblique_sample = 100 * (results_sample(n).combine_fisher_oblique_cardinal_shuffle_sample - ...
                results_sample(n).combine_fisher_oblique_oblique_shuffle_sample) ./ results_sample(n).combine_fisher_oblique_cardinal_shuffle_sample;

        results_sample(n).diff_shuffle_percent_cardinal_median = median(results_sample(n).combine_diff_shuffle_percent_cardinal_sample);
        results_sample(n).diff_shuffle_percent_oblique_median  = median(results_sample(n).combine_diff_shuffle_percent_oblique_sample);

        results_sample(n).diff_shuffle_percent_cardinal_CI     = prctile(results_sample(n).combine_diff_shuffle_percent_cardinal_sample, prc); 
        results_sample(n).diff_shuffle_percent_oblique_CI      = prctile(results_sample(n).combine_diff_shuffle_percent_oblique_sample, prc); 


        %%% average two tasks to get a symmetric measurement of switching
        results_sample(n).combine_diff_info_percent_sample = (results_sample(n).combine_diff_info_percent_cardinal_sample + ...
                                                              results_sample(n).combine_diff_info_percent_oblique_sample) / 2;
        results_sample(n).diff_info_percent_median          = median(results_sample(n).combine_diff_info_percent_sample);
        results_sample(n).diff_info_percent_CI              = prctile(results_sample(n).combine_diff_info_percent_sample,prc); 

        %%%%%% median and CI of I_delta
        %%%%% I_CC_shuffle - I_CC, I_OO_shuffle - I_OO, I_CO_shuffle - I_CO, I_OC_shuffle - I_OC
        results_sample(n).combine_delta_cardinal_cardinal_sample = results_sample(n).combine_fisher_cardinal_cardinal_shuffle_sample - ...
                                                                    results_sample(n).combine_fisher_cardinal_cardinal_sample;  
        results_sample(n).combine_delta_oblique_cardinal_sample  = results_sample(n).combine_fisher_oblique_cardinal_shuffle_sample - ...
                                                                    results_sample(n).combine_fisher_oblique_cardinal_sample;
        results_sample(n).combine_delta_oblique_oblique_sample   = results_sample(n).combine_fisher_oblique_oblique_shuffle_sample - ...
                                                                    results_sample(n).combine_fisher_oblique_oblique_sample;
        results_sample(n).combine_delta_cardinal_oblique_sample  = results_sample(n).combine_fisher_cardinal_oblique_shuffle_sample - ...
                                                                    results_sample(n).combine_fisher_cardinal_oblique_sample;
    
    
        results_sample(n).delta_cardinal_cardinal_median = median(results_sample(n).combine_delta_cardinal_cardinal_sample, 'omitnan');
        results_sample(n).delta_oblique_cardinal_median = median(results_sample(n).combine_delta_oblique_cardinal_sample, 'omitnan');
        results_sample(n).delta_oblique_oblique_median = median(results_sample(n).combine_delta_oblique_oblique_sample, 'omitnan');
        results_sample(n).delta_cardinal_oblique_median = median(results_sample(n).combine_delta_cardinal_oblique_sample, 'omitnan');
        
        results_sample(n).delta_cardinal_cardinal_CI = prctile(results_sample(n).combine_delta_cardinal_cardinal_sample, prc);
        results_sample(n).delta_oblique_cardinal_CI  = prctile(results_sample(n).combine_delta_oblique_cardinal_sample, prc);
        results_sample(n).delta_oblique_oblique_CI   = prctile(results_sample(n).combine_delta_oblique_oblique_sample, prc);
        results_sample(n).delta_cardinal_oblique_CI  = prctile(results_sample(n).combine_delta_cardinal_oblique_sample, prc);
       
        %%%% median and CI of I_delta_pct
        %%%%% 100 * (I_CC_shuffle - I_CC) / I_CC_shuffle 
        %%%%% 100 * (I_OO_shuffle - I_OO) / I_OO_shuffle
        %%%%% 100 * (I_CO_shuffle - I_CO) / I_CO_shuffle
        %%%%% 100 * (I_OC_shuffle - I_OC) / I_OC_shuffle
        results_sample(n).combine_delta_percent_cardinal_cardinal_sample = 100 * (results_sample(n).combine_delta_cardinal_cardinal_sample ./...
                                                                        results_sample(n).combine_fisher_cardinal_cardinal_shuffle_sample);
    
        results_sample(n).combine_delta_percent_oblique_oblique_sample = 100 * (results_sample(n).combine_delta_oblique_oblique_sample ./...
                                                                        results_sample(n).combine_fisher_oblique_oblique_shuffle_sample);

        results_sample(n).combine_delta_percent_cardinal_oblique_sample = 100 * (results_sample(n).combine_delta_cardinal_oblique_sample ./...
                                                                        results_sample(n).combine_fisher_cardinal_oblique_shuffle_sample);

        results_sample(n).combine_delta_percent_oblique_cardinal_sample = 100 * (results_sample(n).combine_delta_oblique_cardinal_sample ./...
                                                                        results_sample(n).combine_fisher_oblique_cardinal_shuffle_sample);

      
        results_sample(n).delta_percent_cardinal_cardinal_median = median(results_sample(n).combine_delta_percent_cardinal_cardinal_sample);
        results_sample(n).delta_percent_oblique_oblique_median = median(results_sample(n).combine_delta_percent_oblique_oblique_sample);
        results_sample(n).delta_percent_cardinal_oblique_median = median(results_sample(n).combine_delta_percent_cardinal_oblique_sample);
        results_sample(n).delta_percent_oblique_cardinal_median = median(results_sample(n).combine_delta_percent_oblique_cardinal_sample);


        results_sample(n).delta_percent_cardinal_cardinal_CI = prctile(results_sample(n).combine_delta_percent_cardinal_cardinal_sample, prc);
        results_sample(n).delta_percent_oblique_oblique_CI = prctile(results_sample(n).combine_delta_percent_oblique_oblique_sample, prc);
        results_sample(n).delta_percent_cardinal_oblique_CI = prctile(results_sample(n).combine_delta_percent_cardinal_oblique_sample, prc);
        results_sample(n).delta_percent_oblique_cardinal_CI = prctile(results_sample(n).combine_delta_percent_oblique_cardinal_sample, prc);

        %%%%% difference in percentage between I_redundancy_within and
        %%%%% I_redundancy_cross
        results_sample(n).combine_diff_delta_percent_cardinal_sample = results_sample(n).combine_delta_percent_cardinal_cardinal_sample - ...
                                                                results_sample(n).combine_delta_percent_cardinal_oblique_sample;

        results_sample(n).combine_diff_delta_percent_oblique_sample = results_sample(n).combine_delta_percent_oblique_oblique_sample - ...
                                                                results_sample(n).combine_delta_percent_oblique_cardinal_sample;

        results_sample(n).diff_delta_percent_cardinal_median    = median(results_sample(n).combine_diff_delta_percent_cardinal_sample);
        results_sample(n).diff_delta_percent_oblique_median     = median(results_sample(n).combine_diff_delta_percent_oblique_sample);

        results_sample(n).diff_delta_percent_cardinal_CI    = prctile(results_sample(n).combine_diff_delta_percent_cardinal_sample, prc);
        results_sample(n).diff_delta_percent_oblique_CI     = prctile(results_sample(n).combine_diff_delta_percent_oblique_sample, prc);
  
    
       
        %%% average two tasks to get a symmetric measurement of switching
        results_sample(n).combine_diff_delta_percent_sample = (results_sample(n).combine_diff_delta_percent_cardinal_sample + ...
                                                              results_sample(n).combine_diff_delta_percent_oblique_sample) / 2;
        results_sample(n).diff_delta_percent_median          = median(results_sample(n).combine_diff_delta_percent_sample);
        results_sample(n).diff_delta_percent_CI              = prctile(results_sample(n).combine_diff_delta_percent_sample,prc); 



        results_sample(n).CI_level = CI_level;
    end
else
        %%%%%% deal with data per coherence %%%%
        for n = 1:nData
        %%%%% median and CI of 4 * 2 fisher information
        results_sample(n).fisher_cardinal_cardinal_median = median(results_sample(n).fisher_cardinal_cardinal_sample, 'omitnan');
        results_sample(n).fisher_oblique_cardinal_median = median(results_sample(n).fisher_oblique_cardinal_sample, 'omitnan');
        results_sample(n).fisher_oblique_oblique_median = median(results_sample(n).fisher_oblique_oblique_sample, 'omitnan');
        results_sample(n).fisher_cardinal_oblique_median = median(results_sample(n).fisher_cardinal_oblique_sample, 'omitnan');
    
        results_sample(n).fisher_cardinal_cardinal_shuffle_median = median(results_sample(n).fisher_cardinal_cardinal_shuffle_sample, 'omitnan');
        results_sample(n).fisher_oblique_cardinal_shuffle_median = median(results_sample(n).fisher_oblique_cardinal_shuffle_sample, 'omitnan');
        results_sample(n).fisher_oblique_oblique_shuffle_median = median(results_sample(n).fisher_oblique_oblique_shuffle_sample, 'omitnan');
        results_sample(n).fisher_cardinal_oblique_shuffle_median = median(results_sample(n).fisher_cardinal_oblique_shuffle_sample, 'omitnan');


        results_sample(n).fisher_cardinal_cardinal_CI = prctile(results_sample(n).fisher_cardinal_cardinal_sample,  prc);
        results_sample(n).fisher_oblique_cardinal_CI = prctile(results_sample(n).fisher_oblique_cardinal_sample,  prc);
        results_sample(n).fisher_oblique_oblique_CI = prctile(results_sample(n).fisher_oblique_oblique_sample,  prc);
        results_sample(n).fisher_cardinal_oblique_CI = prctile(results_sample(n).fisher_cardinal_oblique_sample,  prc);
    
        results_sample(n).fisher_cardinal_cardinal_shuffle_CI = prctile(results_sample(n).fisher_cardinal_cardinal_shuffle_sample,  prc);
        results_sample(n).fisher_oblique_cardinal_shuffle_CI = prctile(results_sample(n).fisher_oblique_cardinal_shuffle_sample,  prc);
        results_sample(n).fisher_oblique_oblique_shuffle_CI = prctile(results_sample(n).fisher_oblique_oblique_shuffle_sample,  prc);
        results_sample(n).fisher_cardinal_oblique_shuffle_CI = prctile(results_sample(n).fisher_cardinal_oblique_shuffle_sample,  prc);
    
    

        %%%% diff between I_cross and I_real as a percentage of I_cross ((I_cross - I_real) / I_cross)
        
        results_sample(n).diff_info_percent_cardinal_sample = 100 * (results_sample(n).fisher_cardinal_oblique_sample - ...
                results_sample(n).fisher_cardinal_cardinal_sample) ./ results_sample(n).fisher_cardinal_oblique_sample;
        results_sample(n).diff_info_percent_oblique_sample = 100 * (results_sample(n).fisher_oblique_cardinal_sample - ...
                results_sample(n).fisher_oblique_oblique_sample) ./ results_sample(n).fisher_oblique_cardinal_sample;

        results_sample(n).diff_info_percent_cardinal_median = median(results_sample(n).diff_info_percent_cardinal_sample);
        results_sample(n).diff_info_percent_oblique_median  = median(results_sample(n).diff_info_percent_oblique_sample);

        results_sample(n).diff_info_percent_cardinal_CI     = prctile(results_sample(n).diff_info_percent_cardinal_sample, prc); 
        results_sample(n).diff_info_percent_oblique_CI      = prctile(results_sample(n).diff_info_percent_oblique_sample, prc); 

        

        %%%% diff between I_cross,shuffle and I_real,shuffle as a percentage of I_cross,shuffle ((I_cross - I_real) / I_cross)
        
        results_sample(n).diff_shuffle_percent_cardinal_sample = 100 * (results_sample(n).fisher_cardinal_oblique_shuffle_sample - ...
                results_sample(n).fisher_cardinal_cardinal_shuffle_sample) ./ results_sample(n).fisher_cardinal_oblique_shuffle_sample;
        results_sample(n).diff_shuffle_percent_oblique_sample = 100 * (results_sample(n).fisher_oblique_cardinal_shuffle_sample - ...
                results_sample(n).fisher_oblique_oblique_shuffle_sample) ./ results_sample(n).fisher_oblique_cardinal_shuffle_sample;

        results_sample(n).diff_shuffle_percent_cardinal_median = median(results_sample(n).diff_shuffle_percent_cardinal_sample);
        results_sample(n).diff_shuffle_percent_oblique_median  = median(results_sample(n).diff_shuffle_percent_oblique_sample);

        results_sample(n).diff_shuffle_percent_cardinal_CI     = prctile(results_sample(n).diff_shuffle_percent_cardinal_sample, prc); 
        results_sample(n).diff_shuffle_percent_oblique_CI      = prctile(results_sample(n).diff_shuffle_percent_oblique_sample, prc); 

        


        %%% average two tasks to get a symmetric measurement of switching
        results_sample(n).diff_info_percent_sample = (results_sample(n).diff_info_percent_cardinal_sample + ...
                                                              results_sample(n).diff_info_percent_oblique_sample) / 2;
        results_sample(n).diff_info_percent_median          = median(results_sample(n).diff_info_percent_sample);
        results_sample(n).diff_info_percent_CI              = prctile(results_sample(n).diff_info_percent_sample,prc); 

        %%%%%% median and CI of I_delta
        %%%%% I_CC_shuffle - I_CC, I_OO_shuffle - I_OO, I_CO_shuffle - I_CO, I_OC_shuffle - I_OC
        results_sample(n).delta_cardinal_cardinal_sample = results_sample(n).fisher_cardinal_cardinal_shuffle_sample - ...
                                                                    results_sample(n).fisher_cardinal_cardinal_sample;  
        results_sample(n).delta_oblique_cardinal_sample  = results_sample(n).fisher_oblique_cardinal_shuffle_sample - ...
                                                                    results_sample(n).fisher_oblique_cardinal_sample;
        results_sample(n).delta_oblique_oblique_sample   = results_sample(n).fisher_oblique_oblique_shuffle_sample - ...
                                                                    results_sample(n).fisher_oblique_oblique_sample;
        results_sample(n).delta_cardinal_oblique_sample  = results_sample(n).fisher_cardinal_oblique_shuffle_sample - ...
                                                                    results_sample(n).fisher_cardinal_oblique_sample;
    
    
        results_sample(n).delta_cardinal_cardinal_median = median(results_sample(n).delta_cardinal_cardinal_sample, 'omitnan');
        results_sample(n).delta_oblique_cardinal_median = median(results_sample(n).delta_oblique_cardinal_sample, 'omitnan');
        results_sample(n).delta_oblique_oblique_median = median(results_sample(n).delta_oblique_oblique_sample, 'omitnan');
        results_sample(n).delta_cardinal_oblique_median = median(results_sample(n).delta_cardinal_oblique_sample, 'omitnan');
        
        results_sample(n).delta_cardinal_cardinal_CI = prctile(results_sample(n).delta_cardinal_cardinal_sample, prc);
        results_sample(n).delta_oblique_cardinal_CI  = prctile(results_sample(n).delta_oblique_cardinal_sample, prc);
        results_sample(n).delta_oblique_oblique_CI   = prctile(results_sample(n).delta_oblique_oblique_sample, prc);
        results_sample(n).delta_cardinal_oblique_CI  = prctile(results_sample(n).delta_cardinal_oblique_sample, prc);
     
        
         %%%% median and CI of I_delta_pct
        %%%%% 100 * (I_CC_shuffle - I_CC) / I_CC_shuffle 
        %%%%% 100 * (I_OO_shuffle - I_OO) / I_OO_shuffle
        %%%%% 100 * (I_CO_shuffle - I_CO) / I_CO_shuffle
        %%%%% 100 * (I_OC_shuffle - I_OC) / I_OC_shuffle
        results_sample(n).delta_percent_cardinal_cardinal_sample = 100 * (results_sample(n).delta_cardinal_cardinal_sample ./...
                                                                        results_sample(n).fisher_cardinal_cardinal_shuffle_sample);
    
        results_sample(n).delta_percent_oblique_oblique_sample = 100 * (results_sample(n).delta_oblique_oblique_sample ./...
                                                                        results_sample(n).fisher_oblique_oblique_shuffle_sample);

        results_sample(n).delta_percent_cardinal_oblique_sample = 100 * (results_sample(n).delta_cardinal_oblique_sample ./...
                                                                        results_sample(n).fisher_cardinal_oblique_shuffle_sample);

        results_sample(n).delta_percent_oblique_cardinal_sample = 100 * (results_sample(n).delta_oblique_cardinal_sample ./...
                                                                        results_sample(n).fisher_oblique_cardinal_shuffle_sample);

      
        results_sample(n).delta_percent_cardinal_cardinal_median = median(results_sample(n).delta_percent_cardinal_cardinal_sample);
        results_sample(n).delta_percent_oblique_oblique_median = median(results_sample(n).delta_percent_oblique_oblique_sample);
        results_sample(n).delta_percent_cardinal_oblique_median = median(results_sample(n).delta_percent_cardinal_oblique_sample);
        results_sample(n).delta_percent_oblique_cardinal_median = median(results_sample(n).delta_percent_oblique_cardinal_sample);


        results_sample(n).delta_percent_cardinal_cardinal_CI = prctile(results_sample(n).delta_percent_cardinal_cardinal_sample, prc);
        results_sample(n).delta_percent_oblique_oblique_CI = prctile(results_sample(n).delta_percent_oblique_oblique_sample, prc);
        results_sample(n).delta_percent_cardinal_oblique_CI = prctile(results_sample(n).delta_percent_cardinal_oblique_sample, prc);
        results_sample(n).delta_percent_oblique_cardinal_CI = prctile(results_sample(n).delta_percent_oblique_cardinal_sample, prc);
    
    



           %%%%% difference in percentage between I_redundancy_within and
        %%%%% I_redundancy_cross
        results_sample(n).diff_delta_percent_cardinal_sample = results_sample(n).delta_percent_cardinal_cardinal_sample - ...
                                                                results_sample(n).delta_percent_cardinal_oblique_sample;

        results_sample(n).diff_delta_percent_oblique_sample = results_sample(n).delta_percent_oblique_oblique_sample - ...
                                                                results_sample(n).delta_percent_oblique_cardinal_sample;

        results_sample(n).diff_delta_percent_cardinal_median    = median(results_sample(n).diff_delta_percent_cardinal_sample);
        results_sample(n).diff_delta_percent_oblique_median     = median(results_sample(n).diff_delta_percent_oblique_sample);

        results_sample(n).diff_delta_percent_cardinal_CI    = prctile(results_sample(n).diff_delta_percent_cardinal_sample, prc);
        results_sample(n).diff_delta_percent_oblique_CI     = prctile(results_sample(n).diff_delta_percent_oblique_sample, prc);
  
    

        %%% average two tasks to get a symmetric measurement of switching
        results_sample(n).diff_delta_percent_sample = (results_sample(n).diff_delta_percent_cardinal_median + ...
                                                              results_sample(n).diff_delta_percent_oblique_median) / 2;
        results_sample(n).diff_delta_percent_median          = median(results_sample(n).diff_delta_percent_sample);
        results_sample(n).diff_delta_percent_CI              = prctile(results_sample(n).diff_delta_percent_sample,prc); 

       

        results_sample(n).CI_level = CI_level;
        end

end
end