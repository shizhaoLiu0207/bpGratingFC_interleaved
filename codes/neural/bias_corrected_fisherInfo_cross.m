function results_crossI = bias_corrected_fisherInfo_cross(X_A, Y_A, X_B,Y_B, X_A_zero, X_B_zero, coherence_level)
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
T_zero_A      = size(X_A_zero, 1);
T_zero_B      = size(X_B_zero, 1);

N_A           = size(X_A, 2); 
N_B           = size(X_B, 2);
assert(N_A == N_B, 'Size of population should be equal under two conditions')
N             = N_A;

if (T_zero_A < (N + 4)) | (T_zero_B < (N + 4)) 
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
%%%% covariance matrix of condition A using zero coherence trials
m_zero_A      = mean(X_A_zero, 1);
C_A_zero      = transpose(X_A_zero - m_zero_A) * (X_A_zero - m_zero_A)  / (T_zero_A - 1);
C_A_zero_diag = diag(diag(C_A_zero));
%%%% covariance matrix of condition B using zero coherence trials
m_zero_B      = mean(X_B_zero, 1);
C_B_zero      = transpose(X_B_zero - m_zero_B) * (X_B_zero - m_zero_B)  / (T_zero_B - 1);
C_B_zero_diag = diag(diag(C_B_zero));

% %%%% covariance matrix of condition A
% C_A = (transpose(X_pos_A - m_pos_A) * (X_pos_A - m_pos_A) + transpose(X_neg_A - m_neg_A) * (X_neg_A - m_neg_A)) / (T1_A + T2_A - 2);
% C_A_diag = diag(diag(C_A));
% %%%% covariance matrix of condition B
% C_B = (transpose(X_pos_B - m_pos_B) * (X_pos_B - m_pos_B) + transpose(X_neg_B - m_neg_B) * (X_neg_B - m_neg_B)) / (T1_B + T2_B - 2);
% C_B_diag = diag(diag(C_B));

%%%%% Naive estimate of IAA: fprime from A, covariance from A (zero coherence trials)
I_AA_naive = fprime_A' * (C_A_zero \ fprime_A);

%%%%% Naive estimate of IBB: fprime from B, covariance from B (zero coherence trials)
I_BB_naive = fprime_B' * (C_B_zero \ fprime_B);

%%%%% Naive estimate of IAB: fprime from A, covariance from B (focus is A, swap the cov matrix)
I_AB_naive = fprime_A' * (C_B_zero \ fprime_A);

%%%%% Naive estimate of IBA: fprime from B, covariance from A (focus is B, swap the cov matrix)
I_BA_naive = fprime_B' * (C_A_zero \ fprime_B);


%%%% Naive estimate of IAA shuffled: fprim from A, diagnol covariance
%%%% matrix from A (zero coherence trials)
I_AA_shuffle_naive = fprime_A' * (C_A_zero_diag \ fprime_A);

%%%% Naive estimate of IBB shuffled: fprim from B, diagnol covariance
%%%% matrix from B (zero coherence trials)
I_BB_shuffle_naive = fprime_B' * (C_B_zero_diag \ fprime_B);

%%%% Naive estimate of IAB shuffled: fprim from A, diagnol covariance matrix from B
I_AB_shuffle_naive = fprime_A' * (C_B_zero_diag \ fprime_A);

%%%% Naive estimate of IBA shuffled: fprim from B, diagnol covariance matrix from A
I_BA_shuffle_naive = fprime_B' * (C_A_zero_diag \ fprime_B);


%%%% bias correction

I_AA_bc = I_AA_naive * (T_zero_A - N - 2) / (T_zero_A - 1) -  (T1_A + T2_A) * N / (T1_A * T2_A * ds^2);
I_BB_bc = I_BB_naive * (T_zero_B - N - 2) / (T_zero_B - 1) -  (T1_B + T2_B) * N / (T1_B * T2_B * ds^2);

I_AB_bc = I_AB_naive * (T_zero_B - N - 2) / (T_zero_B - 1) -  (T1_A + T2_A) * N / (T1_A * T2_A * ds^2);
I_BA_bc = I_BA_naive * (T_zero_A - N - 2) / (T_zero_A - 1) -  (T1_B + T2_B) * N / (T1_B * T2_B * ds^2);

I_AA_shuffle_bc = I_AA_shuffle_naive * (T_zero_A -  3) / (T_zero_A - 1) -  (T1_A + T2_A) * N / (T1_A * T2_A * ds^2);
I_BB_shuffle_bc = I_BB_shuffle_naive * (T_zero_B -  3) / (T_zero_B - 1) -  (T1_B + T2_B) * N / (T1_B * T2_B * ds^2);

I_AB_shuffle_bc = I_AB_shuffle_naive * (T_zero_B -  3) / (T_zero_B - 1) -  (T1_A + T2_A) * N / (T1_A * T2_A * ds^2);
I_BA_shuffle_bc = I_BA_shuffle_naive * (T_zero_A -  3) / (T_zero_A - 1) -  (T1_B + T2_B) * N / (T1_B * T2_B * ds^2);



%%%% variance estimate of I_AA_bc %%%
v = T_zero_A - 1;
p = N;
var_IAA_bc = estimate_cross_variance(I_AA_bc,v,p,N,ds,T1_A,T2_A);

%%%% variance estimate of I_BB_bc %%%
v = T_zero_B - 1;
p = N;
var_IBB_bc = estimate_cross_variance(I_BB_bc,v,p,N,ds,T1_B,T2_B);



%%%% variance estimate of I_AB_bc %%%
v = T_zero_B - 1;
p = N;
var_IAB_bc = estimate_cross_variance(I_AB_bc,v,p,N,ds,T1_A,T2_A);

%%%% variance estimate of I_BA_bc %%%
v = T_zero_A - 1;
p = N;
var_IBA_bc = estimate_cross_variance(I_BA_bc,v,p,N,ds,T1_B,T2_B);


%%%%%%%% variance of shuffled I_AA
v = T_zero_A - 1;
p = 1;
I_AA_unit_naive = fprime_A .^ 2 ./ diag(C_A_zero);
var_IAA_unit      = estimate_cross_variance(I_AA_unit_naive,v,p,N,ds,T1_A,T2_A);
var_IAA_bc_shuffle   = sum(var_IAA_unit);

%%%%%%%% variance of shuffled I_BB
v = T_zero_B - 1;
p = 1;
I_BB_unit_naive = fprime_B .^ 2 ./ diag(C_B_zero);
var_IBB_unit      = estimate_cross_variance(I_BB_unit_naive,v,p,N,ds,T1_B,T2_B);
var_IBB_bc_shuffle   = sum(var_IBB_unit);


%%%%%%%% variance of shuffled I_AB
v = T_zero_B - 1;
p = 1;
I_AB_unit_naive = fprime_A .^ 2 ./ diag(C_B_zero);
var_IAB_unit      = estimate_cross_variance(I_AB_unit_naive,v,p,N,ds,T1_A,T2_A);
var_IAB_bc_shuffle   = sum(var_IAB_unit);

%%%%%%%% variance of shuffled I_BA
v = T_zero_A - 1;
p = 1;
I_BA_unit_naive = fprime_B .^ 2 ./ diag(C_A_zero);
var_IBA_unit      = estimate_cross_variance(I_BA_unit_naive,v,p,N,ds,T1_B,T2_B);
var_IBA_bc_shuffle   = sum(var_IBA_unit);


%%%%%%%% organize result
results_crossI.I_AA_naive = I_AA_naive;
results_crossI.I_BB_naive = I_BB_naive;
results_crossI.I_AB_naive = I_AB_naive;
results_crossI.I_BA_naive = I_BA_naive;

results_crossI.I_AA_bc = I_AA_bc;
results_crossI.I_BB_bc = I_BB_bc;
results_crossI.I_AB_bc = I_AB_bc;
results_crossI.I_BA_bc = I_BA_bc;

results_crossI.I_AA_shuffle_naive = I_AA_shuffle_naive;
results_crossI.I_BB_shuffle_naive = I_BB_shuffle_naive;
results_crossI.I_AB_shuffle_naive = I_AB_shuffle_naive;
results_crossI.I_BA_shuffle_naive = I_BA_shuffle_naive;

results_crossI.I_AA_shuffle_bc = I_AA_shuffle_bc;
results_crossI.I_BB_shuffle_bc = I_BB_shuffle_bc;
results_crossI.I_AB_shuffle_bc = I_AB_shuffle_bc;
results_crossI.I_BA_shuffle_bc = I_BA_shuffle_bc;

results_crossI.var_IAA_bc = var_IAA_bc;
results_crossI.var_IBB_bc = var_IBB_bc;
results_crossI.var_IAB_bc = var_IAB_bc;
results_crossI.var_IBA_bc = var_IBA_bc;

results_crossI.var_IAA_bc_shuffle = var_IAA_bc_shuffle;
results_crossI.var_IBB_bc_shuffle = var_IBB_bc_shuffle;
results_crossI.var_IAB_bc_shuffle = var_IAB_bc_shuffle;
results_crossI.var_IBA_bc_shuffle = var_IBA_bc_shuffle;

function  var_Ibc = estimate_cross_variance(I,v,p,N,ds,T1,T2)

alpha   = 2 / ((v - p) * (v - p - 3));
beta    = (v - p - 1) / ((v - p) * (v - p - 3));
gamma   = (T1 + T2) / (T1 * T2 * ds^2);

var_Ibc = (alpha + 2 * beta) * I .^ 2 + ...
           (6 * alpha + 12 * beta + 4) * gamma * I  + ...
           (3 * alpha + 6 * beta + 2) * gamma ^ 2 * N;

end

end