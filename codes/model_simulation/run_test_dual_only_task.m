clear all
clc
close all
%%
image_task          = 'cardinal';
nTrial              = 100; 
stimulus_contrast   = [0,4];
%%%%% use these priors for dual-task swithcing mode
prior_task_cardinal_list = [0.1:0.2:1];
prior_task_oblique_list  = [0.1:0.2:1];
nCardinal                = numel(prior_task_cardinal_list);
nOblique                 = numel(prior_task_oblique_list) ;

[accuracy_dual, accuracy_single]      = deal(nan*ones(nCardinal, nOblique));  
[accuracy_dual_sem, accuracy_single_sem]      = deal(nan*ones(nCardinal, nOblique));  
for n = 1:nCardinal
    for m = 1:nOblique
        prior_task_cardinal = prior_task_cardinal_list(n);
        prior_task_oblique  = prior_task_oblique_list(m);
        P   = S_Exp_Para('test-dualTask', 'I.stimulus_contrast',stimulus_contrast,'I.image_task',image_task,...
                            'G.prior_task_cardinal', prior_task_cardinal, 'G.prior_task_oblique', prior_task_oblique,...
                            'S.number_repetitions', nTrial);
        dat                     = S_Experiment(P);
        decision                = (dat.O(:,3,end) > 0.5) + 1;
        accuracy_dual(n,m)      = sum(decision == 2) / numel(decision);
        accuracy_dual_sem(n,m)  = sqrt(accuracy_dual(n,m) * (1 - accuracy_dual(n,m)) / numel(decision));
        if abs(prior_task_cardinal + prior_task_oblique - 1) < eps
            
           P   = S_Exp_Para('test-interleaved', 'I.stimulus_contrast',stimulus_contrast,'I.image_task',image_task,...
                            'G.prior_task', [prior_task_cardinal, prior_task_oblique],'S.number_repetitions', nTrial);
            dat                         = S_Experiment(P);
            decision                    = (dat.O(:,3,end) > 0.5) + 1;
            accuracy_single(n,m)        = sum(decision == 2) / numel(decision);
            accuracy_single_sem(n,m)    = sqrt(accuracy_single(n,m) * (1 - accuracy_single(n,m)) / numel(decision));
        end
    end
end
%%
figure;
idx_single = find(~isnan(accuracy_single));
[r,c] = ind2sub(size(accuracy_single),idx_single);
for i = 1:numel(r)
    accuracy_dual_part(i) = accuracy_dual(r(i),c(i));
    accuracy_dual_part_sem(i) = accuracy_dual_sem(r(i),c(i));
end
errorbar(prior_task_cardinal_list(r), accuracy_single(idx_single),accuracy_single_sem(idx_single),'LineWidth',2);hold on
errorbar(prior_task_cardinal_list(r), accuracy_dual_part, accuracy_dual_part_sem,'LineWidth',2)
legend('Single','Dual')
ylabel('Accuracy');
xlabel('Prior_{cardinal}')
set(gca,'fontsize',16)
title('Accuracy for cardinal task, contrast = 4')

figure;
for i = 1:5
    errorbar(prior_task_cardinal_list, accuracy_dual(:,i),accuracy_dual_sem(:,i),'LineWidth',2); hold on
    lgdstr{i} = sprintf('Prior_{oblique} = %.1f',prior_task_oblique_list(i));
end

ylabel('Accuracy');
xlabel('Prior_{cardinal}')
legend(lgdstr)
set(gca,'fontsize',16)
title('Accuracy for cardinal task, contrast = 4 (Dual mode only)')