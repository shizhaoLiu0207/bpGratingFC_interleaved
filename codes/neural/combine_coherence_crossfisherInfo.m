function [dat_fisher_combine,options] = combine_coherence_crossfisherInfo(dat_fisher_cross,options)

if ~isfield(options, 'cohr_limit'), options.cohr_limit = 20; end

taskType_list       = unique({dat_fisher_cross(:).task});
sessionType_list    = unique({dat_fisher_cross(:).sessionType});
timeWinIndex_list   = unique(cell2mat({dat_fisher_cross(:).timeWinIndex}));
sessionStr_list     = unique({dat_fisher_cross(:).sessionStr});

dat_fisher_combine         = struct();
k                   = 1;
nSession            = numel(sessionStr_list);


for i_task = 1:numel(taskType_list)
    for i_type = 1:numel(sessionType_list)
        for i_time = timeWinIndex_list
            taskType        = taskType_list{i_task};
            sessionType     = sessionType_list{i_type};
            
            % linearSlope     = nan * ones(nSession,1);

            for n = 1:nSession
              
                sessionStr = sessionStr_list{n};


                idx             = strcmp({dat_fisher_cross(:).sessionStr}, sessionStr) & ...
                    cell2mat({dat_fisher_cross(:).timeWinIndex}) == i_time & ...
                    strcmp({dat_fisher_cross(:).task},taskType)  & ...
                    strcmp({dat_fisher_cross(:).sessionType},sessionType) & ...
                    cell2mat({dat_fisher_cross(:).coherence_level}) <= options.cohr_limit;
                if sum(idx) == 0
                    continue
                end
                timeWin         = dat_fisher_cross(find(idx,1)).timeWin;
                %%%%%%% combine fisher_cardinal_cardinal information across coherences %%%%%%
                fisher_cardinal_cardinal_bc     = cell2mat({dat_fisher_cross(idx).fisher_cardinal_cardinal_bc});
                var_cardinal_cardial_bc         = cell2mat({dat_fisher_cross(idx).var_cardinal_cardinal_bc});
               
                idx_nan = isnan(fisher_cardinal_cardinal_bc);
                fisher_cardinal_cardinal_bc(idx_nan) = [];
                var_cardinal_cardial_bc(idx_nan) = [];

                w = 1 ./ var_cardinal_cardial_bc / (sum( 1./ var_cardinal_cardial_bc));
                combine_fisher_cardinal_cardinal    =  sum(w .* fisher_cardinal_cardinal_bc);
                combine_var_cardinal_cardinal       =  1 / (sum( 1./ var_cardinal_cardial_bc));
                %%%%%%% combine fisher_oblique_cardinal information across coherences %%%%%%
                fisher_oblique_cardinal_bc     = cell2mat({dat_fisher_cross(idx).fisher_oblique_cardinal_bc});
                var_oblique_cardinal_bc         = cell2mat({dat_fisher_cross(idx).var_oblique_cardinal_bc});

                idx_nan = isnan(fisher_oblique_cardinal_bc);
                fisher_oblique_cardinal_bc(idx_nan) = [];
                var_oblique_cardinal_bc(idx_nan) = [];

                w = 1 ./ var_oblique_cardinal_bc / (sum( 1./ var_oblique_cardinal_bc));
                combine_fisher_oblique_cardinal     =  sum(w .* fisher_oblique_cardinal_bc);
                combine_var_oblique_cardinal        =  1 / (sum( 1./ var_oblique_cardinal_bc));
                %%%%%%% combine fisher_oblique_oblique information across coherences %%%%%%
                fisher_oblique_oblique_bc     = cell2mat({dat_fisher_cross(idx).fisher_oblique_oblique_bc});
                var_oblique_oblique_bc         = cell2mat({dat_fisher_cross(idx).var_oblique_oblique_bc});

                idx_nan = isnan(fisher_oblique_oblique_bc);
                fisher_oblique_oblique_bc(idx_nan) = [];
                var_oblique_oblique_bc(idx_nan) = [];

                w = 1 ./ var_oblique_oblique_bc / (sum( 1./ var_oblique_oblique_bc));
                combine_fisher_oblique_oblique     =  sum(w .* fisher_oblique_oblique_bc);
                combine_var_oblique_oblique        =  1 / (sum( 1./ var_oblique_oblique_bc));
                %%%%%%% combine fisher_cardinal_oblique information across coherences %%%%%%
                fisher_cardinal_oblique_bc     = cell2mat({dat_fisher_cross(idx).fisher_cardinal_oblique_bc});
                var_cardinal_oblique_bc         = cell2mat({dat_fisher_cross(idx).var_cardinal_oblique_bc});
                
                idx_nan = isnan(fisher_cardinal_oblique_bc);
                fisher_cardinal_oblique_bc(idx_nan) = [];
                var_cardinal_oblique_bc(idx_nan) = [];

                w = 1 ./ var_cardinal_oblique_bc / (sum( 1./ var_cardinal_oblique_bc));
                combine_fisher_cardinal_oblique    =  sum(w .* fisher_cardinal_oblique_bc);
                combine_var_cardinal_oblique       =  1 / (sum( 1./ var_cardinal_oblique_bc));
                 %%%%%%% combine fisher_cardinal_cardinal_shuffle information across coherences %%%%%%
                fisher_cardinal_cardinal_shuffle_bc     = cell2mat({dat_fisher_cross(idx).fisher_cardinal_cardinal_shuffle_bc});
                var_cardinal_cardial_shuffle_bc         = cell2mat({dat_fisher_cross(idx).var_cardinal_cardinal_shuffle_bc});

                idx_nan = isnan(fisher_cardinal_cardinal_shuffle_bc);
                fisher_cardinal_cardinal_shuffle_bc(idx_nan) = [];
                var_cardinal_cardial_shuffle_bc(idx_nan) = [];

                w = 1 ./ var_cardinal_cardial_shuffle_bc / (sum( 1./ var_cardinal_cardial_shuffle_bc));
                combine_fisher_cardinal_cardinal_shuffle    =  sum(w .* fisher_cardinal_cardinal_shuffle_bc);
                combine_var_cardinal_cardinal_shuffle       =  1 / (sum( 1./ var_cardinal_cardial_shuffle_bc));
                %%%%%%% combine fisher_oblique_cardinal_shuffle information across coherences %%%%%%
                fisher_oblique_cardinal_shuffle_bc     = cell2mat({dat_fisher_cross(idx).fisher_oblique_cardinal_shuffle_bc});
                var_oblique_cardinal_shuffle_bc         = cell2mat({dat_fisher_cross(idx).var_oblique_cardinal_shuffle_bc});

                idx_nan = isnan(fisher_oblique_cardinal_shuffle_bc);
                fisher_oblique_cardinal_shuffle_bc(idx_nan) = [];
                var_oblique_cardinal_shuffle_bc(idx_nan) = [];

                w = 1 ./ var_oblique_cardinal_shuffle_bc / (sum( 1./ var_oblique_cardinal_shuffle_bc));
                combine_fisher_oblique_cardinal_shuffle     =  sum(w .* fisher_oblique_cardinal_shuffle_bc);
                combine_var_oblique_cardinal_shuffle        =  1 / (sum( 1./ var_oblique_cardinal_shuffle_bc));
                %%%%%%% combine fisher_oblique_oblique_shuffle information across coherences %%%%%%
                fisher_oblique_oblique_shuffle_bc     = cell2mat({dat_fisher_cross(idx).fisher_oblique_oblique_shuffle_bc});
                var_oblique_oblique_shuffle_bc         = cell2mat({dat_fisher_cross(idx).var_oblique_oblique_shuffle_bc});

                idx_nan = isnan(fisher_oblique_oblique_shuffle_bc);
                fisher_oblique_oblique_shuffle_bc(idx_nan) = [];
                var_oblique_oblique_shuffle_bc(idx_nan) = [];


                w = 1 ./ var_oblique_oblique_shuffle_bc / (sum( 1./ var_oblique_oblique_shuffle_bc));
                combine_fisher_oblique_oblique_shuffle     =  sum(w .* fisher_oblique_oblique_shuffle_bc);
                combine_var_oblique_oblique_shuffle        =  1 / (sum( 1./ var_oblique_oblique_shuffle_bc));
                %%%%%%% combine fisher_cardinal_oblique_shuffle information across coherences %%%%%%
                fisher_cardinal_oblique_shuffle_bc     = cell2mat({dat_fisher_cross(idx).fisher_cardinal_oblique_shuffle_bc});
                var_cardinal_oblique_shuffle_bc         = cell2mat({dat_fisher_cross(idx).var_cardinal_oblique_shuffle_bc});

                idx_nan = isnan(fisher_cardinal_oblique_shuffle_bc);
                fisher_cardinal_oblique_shuffle_bc(idx_nan) = [];
                var_cardinal_oblique_shuffle_bc(idx_nan) = [];

                w = 1 ./ var_cardinal_oblique_shuffle_bc / (sum( 1./ var_cardinal_oblique_shuffle_bc));
                combine_fisher_cardinal_oblique_shuffle    =  sum(w .* fisher_cardinal_oblique_shuffle_bc);
                combine_var_cardinal_oblique_shuffle       =  1 / (sum( 1./ var_cardinal_oblique_shuffle_bc));

              
                




                dat_fisher_combine(k).sessionStr                    = sessionStr;
                dat_fisher_combine(k).taskType                      = taskType;
                dat_fisher_combine(k).sessionType                   = sessionType;
                dat_fisher_combine(k).timeWinIndex                  = i_time;
                dat_fisher_combine(k).timeWin                       = timeWin;
                dat_fisher_combine(k).nUnit                         = dat_fisher_cross(1).nNeuron;
                dat_fisher_combine(k).combine_fisher_cardinal_cardinal           = combine_fisher_cardinal_cardinal;
                dat_fisher_combine(k).combine_fisher_cardinal_cardinal_std       = sqrt(combine_var_cardinal_cardinal);

                dat_fisher_combine(k).combine_fisher_oblique_cardinal           = combine_fisher_oblique_cardinal;
                dat_fisher_combine(k).combine_fisher_oblique_cardinal_std       = sqrt(combine_var_oblique_cardinal);

                dat_fisher_combine(k).combine_fisher_oblique_oblique            = combine_fisher_oblique_oblique;
                dat_fisher_combine(k).combine_fisher_oblique_oblique_std       = sqrt(combine_var_oblique_oblique);

                dat_fisher_combine(k).combine_fisher_cardinal_oblique           = combine_fisher_cardinal_oblique;
                dat_fisher_combine(k).combine_fisher_cardinal_oblique_std       = sqrt(combine_var_cardinal_oblique);


               
                dat_fisher_combine(k).combine_fisher_cardinal_cardinal_shuffle           = combine_fisher_cardinal_cardinal_shuffle;
                dat_fisher_combine(k).combine_fisher_cardinal_cardinal_shuffle_std       = sqrt(combine_var_cardinal_cardinal_shuffle);

                dat_fisher_combine(k).combine_fisher_oblique_cardinal_shuffle           = combine_fisher_oblique_cardinal_shuffle;
                dat_fisher_combine(k).combine_fisher_oblique_cardinal_shuffle_std       = sqrt(combine_var_oblique_cardinal_shuffle);

                dat_fisher_combine(k).combine_fisher_oblique_oblique_shuffle            = combine_fisher_oblique_oblique_shuffle;
                dat_fisher_combine(k).combine_fisher_oblique_oblique_shuffle_std       =  sqrt(combine_var_oblique_oblique_shuffle);

                dat_fisher_combine(k).combine_fisher_cardinal_oblique_shuffle           = combine_fisher_cardinal_oblique_shuffle;
                dat_fisher_combine(k).combine_fisher_cardinal_oblique_shuffle_std       = sqrt(combine_var_cardinal_oblique_shuffle);

            
                k = k+1;
            end
        end
    end
end