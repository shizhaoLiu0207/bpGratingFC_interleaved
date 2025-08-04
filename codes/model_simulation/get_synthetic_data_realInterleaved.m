function  data_out = get_synthetic_data_realInterleaved(synthData_use, timeBinIndex)

nRepeat                 = size(synthData_use(1).X_response,1);
idx_data                = find(cell2mat({synthData_use(:).timebinIndex}) == timeBinIndex); % find data correspond to this timebin
nCondition              = numel(idx_data); % number of stimulus conditions


data_out                = struct();
% tmp_original            = {synthData_use(idx_data).X_response_original};
% tmp_shift               = {synthData_use(idx_data).X_response_shift};
% data_out.spikeCount     = [cat(1,tmp_original{:}); cat(1,tmp_shift{:})];
tmp                     = {synthData_use(:).X_response};
data_out.spikeCount     = cat(1, tmp{:});


nTrial                  = size(data_out.spikeCount,1);
data_out.signal         = zeros(nTrial,1);
data_out.orientation    = zeros(nTrial,1);
data_out.behavCorrect   = zeros(nTrial,1);
data_out.timebin        = synthData_use(idx_data(1)).timebin;



for i = 1:nCondition
    idx_trial = [(i-1) * nRepeat  + 1 : i * nRepeat];
    contrast_signed = synthData_use(idx_data(i)).contrast_signed;
    image_task      =  synthData_use(idx_data(i)).image_task;


    data_out.signal(idx_trial) = contrast_signed;


    switch image_task
        case 'cardinal'
            if contrast_signed > 0
                % 90 degree
               data_out.orientation(idx_trial) = 90;
            elseif contrast_signed < 0
                % 0 degree
                data_out.orientation(idx_trial) = 0;
            elseif contrast_signed  == 0
                % randomly assign
                data_out.orientation(idx_trial) = randsample([0,90], numel(idx_trial), true);
            end
        case 'oblique'      
            if contrast_signed > 0
                % 135 degree
                data_out.orientation(idx_trial) = 135;
            elseif contrast_signed < 0
                % 45 degree
                data_out.orientation(idx_trial) = 45;
            elseif contrast_signed  == 0
                % randomly assign
                data_out.orientation(idx_trial) = randsample([45,135], numel(idx_trial), true);
            end
    end

    %%%% add behavior

    switch image_task
        case 'cardinal'
            ori_list = [0,90];
        case 'oblique'
            ori_list = [45,135];
    end
    choiceOri = ori_list(synthData_use(idx_data(i)).decision);
    data_out.behavCorrect(idx_trial) = choiceOri' == data_out.orientation(idx_trial);
end


    
% end
% for i  = 1:nCondition
% 
%     idx_trial =  [(i-1) * nRepeat  + 1 : i * nRepeat] + nCondition * nRepeat;
%     contrast = synthData_use(idx_data(i)).stimulus_contrast;
%     data_out.signal(idx_trial) = contrast(1) - contrast(2);
%     if contrast(1) >= contrast(2)
%         data_out.orientation(idx_trial) = 45;
%     else
%         data_out.orientation(idx_trial) = 135;
%     end
% end
% 







end