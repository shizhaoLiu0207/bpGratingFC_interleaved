clear all
clc
close all
%%
saveName = '../../results/neural/simulation_results_estimateCrossFisher';
doRun = 0;
if doRun
    N_list   = [2, 10, 30];
    rho_list = [0.6:0.2:1];
    T_list   = [50,100,200,500,1000,2000];

    k = 1;
    for N = N_list
        N_simulation = 1000;


            ds = 5;
            mu1_cardinal = randsample([10:30],N,'true')';
            mu2_cardinal = randsample([20:40],N,'true')';
            mu0_cardinal = (mu1_cardinal + mu2_cardinal) / 2;

            mu1_oblique  = randsample([20:40],N,'true')';
            mu2_oblique  = randsample([10:30],N,'true')';
            mu0_oblique  = (mu1_oblique + mu2_oblique) / 2;

            %mu2 = mu1+ randn(size(mu1));
            
      
            sigma_1_cardinal = sqrt(mu1_cardinal);
            sigma_2_cardinal = sqrt(mu2_cardinal);
            sigma_0_cardinal = sqrt(mu0_cardinal);
            % sigma_cardinal = sigma_1_cardinal + sigma_2_cardinal;

            sigma_1_oblique = sqrt(mu1_oblique);
            sigma_2_oblique = sqrt(mu2_oblique);
            sigma_0_oblique = sqrt(mu0_oblique);
            % sigma_oblique = sigma_1_oblique + sigma_2_oblique;

            [C_cardinal, C_diag_cardinal] = generate_C(sigma_0_cardinal, N);
            [C_oblique, C_diag_oblique]  = generate_C(sigma_0_oblique, N);
         


            f_prime_cardinal  = (mu2_cardinal - mu1_cardinal) / (2 * ds);
            f_prime_oblique   = (mu2_oblique - mu1_oblique) / (2 * ds);
            
            trueI_cardinal_oblique = f_prime_cardinal' * (C_oblique \ f_prime_cardinal);
            trueI_oblique_cardinal =  f_prime_oblique' * (C_cardinal \ f_prime_oblique);

            trueI_cardinal_oblique_shuffle = f_prime_cardinal' * (C_diag_oblique \ f_prime_cardinal);
            trueI_oblique_cardinal_shuffle =  f_prime_oblique' * (C_diag_cardinal \ f_prime_oblique);
            % trueI = f_prime' * (C \ f_prime);
            % trueI_shuffle = f_prime' * (C_diag \ f_prime);



        for rho = rho_list

            %%%%%
            N_T = numel(T_list);
            %[I_real_naive, I_real_bc, I_real_bc_new,I_shuffle_naive, I_shuffle_bc, I_shuffle_bc_new,var_I_real, var_I_shuffle] = deal(nan * ones(N_simulation, N_T));
            % [var_I_real_paper,var_I_shuffle_paper] = deal(nan * ones(N_simulation, N_T));
            [I_real_cardinal_oblique_naive, I_real_cardinal_oblique_bc] = deal(nan * ones(N_simulation, N_T));
            [I_shuffle_cardinal_oblique_naive, I_shuffle_cardinal_oblique_bc] = deal(nan * ones(N_simulation, N_T));
            [I_real_oblique_cardinal_naive, I_real_oblique_cardinal_bc] = deal(nan * ones(N_simulation, N_T));
            [I_shuffle_oblique_cardinal_naive, I_shuffle_oblique_cardinal_bc] = deal(nan * ones(N_simulation, N_T));
            for t = 1:numel(T_list)


                T1 = T_list(t);
                T2 = T_list(t) * rho;
                T0 = round(T1 / 2);

                for n = 1:N_simulation
                    
                    X1_cardinal = mvnrnd(mu1_cardinal, C_cardinal, T1);
                    X2_cardinal = mvnrnd(mu2_cardinal, C_cardinal, T2);
                    X0_cardinal = mvnrnd(mu0_cardinal, C_cardinal, T0);


                    X1_oblique = mvnrnd(mu1_oblique, C_oblique, T1);
                    X2_oblique = mvnrnd(mu2_oblique, C_oblique, T2);
                    X0_oblique = mvnrnd(mu0_oblique, C_oblique, T0);

                    X_cardinal = [X1_cardinal;X2_cardinal];
                    Y_cardinal = [ones(T1,1); -1 * ones(T2,1)];

                    X_oblique = [X1_oblique;X2_oblique];
                    Y_oblique = [ones(T1,1); -1 * ones(T2,1)];

                    results_crossI = bias_corrected_fisherInfo_cross(X_cardinal, Y_cardinal, X_oblique,Y_oblique, X0_cardinal, X0_oblique, ds);

                    

                    if results_crossI.enoughTrials

                        I_real_cardinal_oblique_naive(n,t)  = results_crossI.I_AB_naive;
                        I_real_cardinal_oblique_bc(n,t)     = results_crossI.I_AB_bc;

                        I_real_oblique_cardinal_naive(n,t)  = results_crossI.I_BA_naive;
                        I_real_oblique_cardinal_bc(n,t)     = results_crossI.I_BA_bc;


                        I_shuffle_cardinal_oblique_naive(n,t)   = results_crossI.I_AB_shuffle_naive;
                        I_shuffle_cardinal_oblique_bc(n,t)      = results_crossI.I_AB_shuffle_bc;

                        I_shuffle_oblique_cardinal_naive(n,t)   = results_crossI.I_BA_shuffle_naive;
                        I_shuffle_oblique_cardinal_bc(n,t)      = results_crossI.I_BA_shuffle_bc;


                   
                    end
                end
            end
            results_simulation(k).N                 = N;
            results_simulation(k).T_list            = T_list;
            results_simulation(k).rho               = rho;

            results_simulation(k).I_real_cardinal_oblique_naive         = I_real_cardinal_oblique_naive;
            results_simulation(k).I_real_cardinal_oblique_bc            = I_real_cardinal_oblique_bc;
            results_simulation(k).I_real_oblique_cardinal_naive         = I_real_oblique_cardinal_naive;
            results_simulation(k).I_real_oblique_cardinal_bc            = I_real_oblique_cardinal_bc;

            results_simulation(k).I_shuffle_cardinal_oblique_naive      = I_shuffle_cardinal_oblique_naive;
            results_simulation(k).I_shuffle_cardinal_oblique_bc         = I_shuffle_cardinal_oblique_bc;
            results_simulation(k).I_shuffle_oblique_cardinal_naive      = I_shuffle_oblique_cardinal_naive;
            results_simulation(k).I_shuffle_oblique_cardinal_bc         = I_shuffle_oblique_cardinal_bc;
  
            results_simulation(k).trueI_cardinal_oblique                = trueI_cardinal_oblique;
            results_simulation(k).trueI_oblique_cardinal                = trueI_oblique_cardinal;

            results_simulation(k).trueI_cardinal_oblique_shuffle        = trueI_cardinal_oblique_shuffle;
            results_simulation(k).trueI_oblique_cardinal_shuffle        = trueI_oblique_cardinal_shuffle;
 
            k = k+1;
        end
        save(saveName,'results_simulation')
   
        
    end
else
    load(saveName);
end


%% visualization of I_bc as a function of T under each condition

N_list      = unique(cell2mat({results_simulation(:).N}));
rho_list    = unique(cell2mat({results_simulation(:).rho}));


for i_N = 1:numel(N_list)
    figure
    for i_rho = 1:numel(rho_list)
        subplot(1,numel(rho_list),i_rho)
        idx = cell2mat({results_simulation(:).N}) == N_list(i_N) & ...
                cell2mat({results_simulation(:).rho}) == rho_list(i_rho);
      
        h = plot_I(results_simulation(idx).I_real_cardinal_oblique_naive, ...
                    results_simulation(idx).I_real_cardinal_oblique_bc,...
                    results_simulation(idx).I_shuffle_cardinal_oblique_naive,...
                    results_simulation(idx).I_shuffle_cardinal_oblique_bc,...
                    results_simulation(idx).trueI_cardinal_oblique, ...
                    results_simulation(idx).trueI_cardinal_oblique_shuffle,...
                    results_simulation(idx).T_list);

%         plot_I_new(results_simulation(idx).I_real_bc,results_simulation(idx).I_real_bc_new,...
%             results_simulation(idx).I_shuffle_bc,results_simulation(idx).I_shuffle_bc_new,...
%             results_simulation(idx).trueI, results_simulation(idx).trueI_shuffle,T_list);


        xlabel('T1');
        ylabel('Fisher Information')
        title(sprintf('T2 = %.2f * T1',rho_list(i_rho)));
    end
    sgtitle(sprintf('Cardinal-oblique; N = %d',N_list(i_N)),'fontsize', 25)
end


N_list      = unique(cell2mat({results_simulation(:).N}));
rho_list    = unique(cell2mat({results_simulation(:).rho}));


for i_N = 1:numel(N_list)
    figure
    for i_rho = 1:numel(rho_list)
        subplot(1,numel(rho_list),i_rho)
        idx = cell2mat({results_simulation(:).N}) == N_list(i_N) & ...
                cell2mat({results_simulation(:).rho}) == rho_list(i_rho);
      
        h = plot_I(results_simulation(idx).I_real_oblique_cardinal_naive, ...
                    results_simulation(idx).I_real_oblique_cardinal_bc,...
                    results_simulation(idx).I_shuffle_oblique_cardinal_naive,...
                    results_simulation(idx).I_shuffle_oblique_cardinal_bc,...
                    results_simulation(idx).trueI_oblique_cardinal, ...
                    results_simulation(idx).trueI_oblique_cardinal_shuffle,...
                    results_simulation(idx).T_list);

%         plot_I_new(results_simulation(idx).I_real_bc,results_simulation(idx).I_real_bc_new,...
%             results_simulation(idx).I_shuffle_bc,results_simulation(idx).I_shuffle_bc_new,...
%             results_simulation(idx).trueI, results_simulation(idx).trueI_shuffle,T_list);


        xlabel('T1');
        ylabel('Fisher Information')
        title(sprintf('T2 = %.2f * T1',rho_list(i_rho)));
    end
    sgtitle(sprintf('Oblique-cardinal; N = %d',N_list(i_N)),'fontsize', 25)
end
%% helper functions
function [C, C_diag] = generate_C(sigma, N)
flag_pos_semi = 0;
while ~flag_pos_semi
    A =  2 * randn(N);
    rho_m = A' * A / N^2;
    rho_m([1:N+1:N^2]) = 1;
    C = sqrt(sigma)' .* rho_m .* sqrt(sigma);
    if all(eig(C) > 0)
        flag_pos_semi = 1;
    end
end
C_diag = diag(diag(C));
end

function h = plot_I(I_real_naive,I_real_bc,I_shuffle_naive,I_shuffle_bc,trueI, trueI_shuffle,T_list)

I_real_naive_avg    = mean(I_real_naive,1);
I_real_bc_avg       = mean(I_real_bc,1);
I_real_naive_sem    = std(I_real_naive ,[], 1) / sqrt(size(I_real_naive,1));
I_real_bc_sem       = std(I_real_bc, [],1) / sqrt(size(I_real_bc,1));

I_shuffle_naive_avg    = mean(I_shuffle_naive,1);
I_shuffle_bc_avg       = mean(I_shuffle_bc,1);
I_shuffle_naive_sem    = std(I_shuffle_naive ,[], 1) / sqrt(size(I_shuffle_naive,1));
I_shuffle_bc_sem       = std(I_shuffle_bc, [],1) / sqrt(size(I_shuffle_bc,1));


x    = log(T_list);
h(1) = plot(x, I_real_naive_avg, 'Linewidth',3); hold on
h(2) = plot(x, I_real_bc_avg , 'LineWidth',3);
h(3) = plot(x,trueI * ones(numel(T_list),1),'linestyle','-','color','black','linewidth',3);
h(4) = plot(x, I_shuffle_naive_avg, 'Linewidth',3, 'linestyle','--','color',h(1).Color);
h(5) = plot(x, I_shuffle_bc_avg, 'Linewidth',3 , 'linestyle', '--','color',h(2).Color);
h(6) = plot(x, trueI_shuffle * ones(numel(T_list),1),'linestyle','--','color','black','linewidth',3);

fill([x, fliplr(x)], [I_real_naive_avg - I_real_naive_sem, fliplr(I_real_naive_avg + I_real_naive_sem)], h(1).Color,...
    'FaceAlpha',0.5,'edgecolor','none')
fill([x, fliplr(x)], [I_real_bc_avg - I_real_bc_sem, fliplr(I_real_bc_avg + I_real_bc_sem)], h(2).Color,...
    'FaceAlpha',0.5,'edgecolor','none')


fill([x, fliplr(x)], [I_shuffle_naive_avg - I_shuffle_naive_sem, fliplr(I_shuffle_naive_avg + I_shuffle_naive_sem)], h(1).Color,...
    'FaceAlpha',0.5,'edgecolor','none')
fill([x, fliplr(x)], [I_shuffle_bc_avg - I_shuffle_bc_sem, fliplr(I_shuffle_bc_avg + I_shuffle_bc_sem)], h(2).Color,...
    'FaceAlpha',0.5,'edgecolor','none')

set(gca,'xtick',log(T_list),'xticklabels',T_list,'fontsize',20);
xlabel('T');
legend('I_real_naive','I_real_bc','true_I_real','I_shuffle_naive','I_shuffle_bc','true_I_shuffle','interpreter','none')
end