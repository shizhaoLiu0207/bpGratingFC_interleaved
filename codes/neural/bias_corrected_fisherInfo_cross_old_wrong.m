function results_crossI = bias_corrected_fisherInfo_cross_old_wrong(X_A, Y_A, X_B,Y_B, coherence_level)
%%%%% Nov-20-2024, Shizhao: This script is using the best way to compute
%%%%% cross information. The most appropriate way is to use covariance from
%%%%% zero coherence trials, which controls stimulus and thus variance.
%%%%% The version using covariance matrix from zero coherence trials is
%%%%% bias_corrected_fisherInfo_cross.m

% Reference Kanitscheider, I., Coen-Cagli, R., Kohn, A., & Pouget, A.
% (2015). Measuring Fisher Information Accurately in Correlated Neural Populations. PLoS
% Computational Biology, 11(6), 1-27. https://doi.org/10.1371/journal.pcbi.1004218
% extend two equal T to different trial numbers for two categories, T1 and
% T2.

% Bias-corrected "cross" information. Using fprime for one dataset but cov
% from the other dataset



% X_A:  neural response of two categories in condition A
% Y_A: label of two categories in condition A
% X_B:  neural response of two categories in condition B
% Y_B: label of two categories in condition B
% coherence level: self-explanatory
T1_A          = sum(Y_A == 1);
T2_A          = sum(Y_A == -1);
T1_B          = sum(Y_B == 1);
T2_B          = sum(Y_B == -1);

N_A           = size(X_A, 2); 
N_B           = size(X_B, 2);
assert(N_A == N_B, 'Size of population should be equal under two conditions')
N             = N_A;

if ((T1_A + T2_A) < N + 5) | (T1_A < (N+1) / 2) | (T2_A < (N+1) / 2) | ...
        ((T1_B + T2_B) < N + 5) | (T1_B < (N+1) / 2) | (T2_B < (N+1) / 2) 
    results_crossI.I_AB_naive = nan;
    results_crossI.I_BA_naive = nan;
    results_crossI.I_AB_bc = nan;
    results_crossI.I_BA_bc = nan;
    results_crossI.var_IAB_bc = nan;
    results_crossI.var_IBA_bc = nan;
    results_crossI.enoughTrials = 0;
    return
end
results_crossI.enoughTrials = 1;

ds          = 2 * coherence_level; 

%%%% fprime of condition A
X_pos_A       = X_A(Y_A == 1,:);
X_neg_A       = X_A(Y_A == -1,:);
m_pos_A       = mean(X_pos_A,1);
m_neg_A       = mean(X_neg_A,1);
fprime_A      = transpose(m_pos_A - m_neg_A) / ds; 
%%%% fprime of condition B
X_pos_B       = X_B(Y_B == 1,:);
X_neg_B       = X_B(Y_B == -1,:);
m_pos_B       = mean(X_pos_B,1);
m_neg_B       = mean(X_neg_B,1);
fprime_B      = transpose(m_pos_B - m_neg_B) / ds; 
%%%% covariance matrix of condition A
C_A = (transpose(X_pos_A - m_pos_A) * (X_pos_A - m_pos_A) + transpose(X_neg_A - m_neg_A) * (X_neg_A - m_neg_A)) / (T1_A + T2_A - 2);
C_A_diag = diag(diag(C_A));
%%%% covariance matrix of condition B
C_B = (transpose(X_pos_B - m_pos_B) * (X_pos_B - m_pos_B) + transpose(X_neg_B - m_neg_B) * (X_neg_B - m_neg_B)) / (T1_B + T2_B - 2);
C_B_diag = diag(diag(C_B));
%%%%% Naive estimate of IAB: fprime from A, covariance from B (focus is A, swap the cov matrix)
I_AB_naive = fprime_A' * (C_B \ fprime_A);

%%%%% Naive estimate of IBA: fprime from B, covariance from A (focus is B, swap the cov matrix)
I_BA_naive = fprime_B' * (C_A \ fprime_B);

%%%% Naive estimate of IAB shuffled: fprim from A, diagnol covariance matrix from B
I_AB_shuffle_naive = fprime_A' * (C_B_diag \ fprime_A);

%%%% Naive estimate of IBA shuffled: fprim from B, diagnol covariance matrix from A
I_BA_shuffle_naive = fprime_B' * (C_A_diag \ fprime_B);

%%%% bias correction
I_AB_bc = I_AB_naive * (T1_B + T2_B - N - 3) / (T1_B + T2_B - 2) -  (T1_A + T2_A) * N / (T1_A * T2_A * ds^2);
I_BA_bc = I_BA_naive * (T1_A + T2_A - N - 3) / (T1_A + T2_A - 2) -  (T1_B + T2_B) * N / (T1_B * T2_B * ds^2);

I_AB_shuffle_bc = I_AB_shuffle_naive * (T1_B + T2_B -  4) / (T1_B + T2_B - 2) -  (T1_A + T2_A) * N / (T1_A * T2_A * ds^2);
I_BA_shuffle_bc = I_BA_shuffle_naive * (T1_A + T2_A -  4) / (T1_A + T2_A - 2) -  (T1_B + T2_B) * N / (T1_B * T2_B * ds^2);


%%%% variance estimate of I_AB_bc %%%
v = T1_B + T2_B - 2;
p = N;
var_IAB_bc = estimate_cross_variance(I_AB_bc,v,p,N,ds,T1_A,T2_A);

%%%% variance estimate of I_BA_bc %%%
v = T1_A + T2_A - 2;
p = N;
var_IBA_bc = estimate_cross_variance(I_BA_bc,v,p,N,ds,T1_B,T2_B);

%%%%% test: use variance and fprime from task A, and correlation matrix
%%%%% from task B
var_A = ( sum((X_pos_A - m_pos_A) .^ 2, 1) +  sum( (X_neg_A - m_neg_A) .^ 2, 1) ) / (T1_A + T2_A - 2);
var_B = ( sum((X_pos_B - m_pos_B) .^ 2, 1) +  sum((X_neg_B - m_neg_B) .^ 2, 1) ) / (T1_B + T2_B - 2);
rho_A = C_A ./ ( transpose(sqrt(var_A)) * sqrt(var_A));
rho_B = C_B ./ ( transpose(sqrt(var_B)) * sqrt(var_B));

C_cross_AB = diag(sqrt(var_A)) *  rho_B * diag(sqrt(var_A));
C_cross_BA = diag(sqrt(var_B)) *  rho_A * diag(sqrt(var_B));

%%%%% Naive estimate of IAB: fprime from A, var from A, correlation from B (focus is A, swap the correlation matrix)
I_AB_naive_corr = fprime_A' * (C_cross_AB \ fprime_A);
%%%%% Naive estimate of IBA: fprime from B, covariance from A (focus is B, swap the cov matrix)
I_BA_naive_corr = fprime_B' * (C_cross_BA \ fprime_B);

I_AB_bc_corr = I_AB_naive_corr * (T1_A + T2_A - N - 3) / (T1_A + T2_A - 2) -  (T1_A + T2_A) * N / (T1_A * T2_A * ds^2);
I_BA_bc_corr = I_BA_naive_corr * (T1_B + T2_B - N - 3) / (T1_B + T2_B - 2) -  (T1_B + T2_B) * N / (T1_B * T2_B * ds^2);

%%%%%
results_crossI.I_AB_naive = I_AB_naive;
results_crossI.I_BA_naive = I_BA_naive;
results_crossI.I_AB_bc = I_AB_bc;
results_crossI.I_BA_bc = I_BA_bc;
results_crossI.I_AB_shuffle_naive = I_AB_shuffle_naive;
results_crossI.I_BA_shuffle_naive = I_BA_shuffle_naive;
results_crossI.I_AB_shuffle_bc = I_AB_shuffle_bc;
results_crossI.I_BA_shuffle_bc = I_BA_shuffle_bc;

results_crossI.var_IAB_bc = var_IAB_bc;
results_crossI.var_IBA_bc = var_IBA_bc;

results_crossI.I_AB_naive_corr = I_AB_naive_corr;
results_crossI.I_BA_naive_corr = I_BA_naive_corr;

results_crossI.I_AB_bc_corr = I_AB_bc_corr;
results_crossI.I_BA_bc_corr = I_BA_bc_corr;

function  var_Ibc = estimate_cross_variance(I,v,p,N,ds,T1,T2)

alpha   = 2 / ((v - p) * (v - p - 3));
beta    = (v - p - 1) / ((v - p) * (v - p - 3));
gamma   = (T1 + T2) / (T1 * T2 * ds^2);

var_Ibc = (alpha + 2 * beta) * I .^ 2 + ...
           (6 * alpha + 12 * beta + 4) * gamma * I  + ...
           (3 * alpha + 6 * beta + 2) * gamma ^ 2 * N;

end

end