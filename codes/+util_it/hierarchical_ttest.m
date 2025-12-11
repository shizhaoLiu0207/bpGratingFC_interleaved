function lmm_results = hierarchical_ttest(rate, group, session)
% Inputs you should have as column vectors (length = number of neuron×trial samples):
% rate         : mean firing rate over your analysis window (Hz)
% group        : 0/1 or char/categorical indicating neuron type (e.g., 'A','B')
% session      : session ID (string/int -> categorical)
% neuron       : neuron ID (optional; only if you can track putative repeats)
% monkey       : optional higher level (if you have multiple animals)

% ---------- Build table ----------
T = table();
T.rate    = rate(:);
%T.lograte = log(rate(:) + 1e-6);  % small epsilon to avoid log(0)
T.group   = categorical(group(:));      % e.g., 'A','B'
T.session = categorical(session(:));
if exist('neuron','var'), T.neuron = categorical(neuron(:)); end
if exist('monkey','var'), T.monkey = categorical(monkey(:)); end

% ---------- Specify the hierarchical structure ----------
% Random intercepts per session; (optional) crossed random intercepts for neuron ID;
% (optional) random intercepts for monkey.
re = "(1|session)";
if ismember('neuron', T.Properties.VariableNames)
    re = re + " + (1|neuron)";   % crossed REs: session and neuron
end
if ismember('monkey', T.Properties.VariableNames)
    re = re + " + (1|monkey)";
end

% Fixed effect for group (this is your “t-test” effect):
formula = "rate ~ 1 + group + " + re;

% ---------- Fit model ----------
lme = fitlme(T, formula, 'DummyVarCoding','effects');  % effect coding → group coeff = A vs B mean diff

% ---------- Test A vs B ----------
% With effects coding and 2 levels, the group coefficient equals half the A–B difference.
% Easier: use coefTest with a contrast on the group term as returned.
disp(lme)
[beta,~,stats] = fixedEffects(lme,'DFMethod','satterthwaite');
ci = coefCI(lme,'DFMethod','satterthwaite');

% Identify which row is the group effect:
fxNames = lme.CoefficientNames;  % {'(Intercept)','group_A'} or similar with effects coding
ix = find(contains(fxNames,'group'));

delta       = beta(ix);              % effect on log-rate
delta_CI    = ci(ix,:);             % CI on log-rate scale
pval        = stats.pValue(ix);
tval        = stats.tStat(ix);
df          = stats.DF(ix); 
% Convert to % difference on the original rate scale (approx: exp(delta)-1):
delta_pct  = delta * 100;
delta_pct_CI = delta_CI * 100;

lmm_results.delta_pct       = delta_pct;
lmm_results.delta_pct_CI    = delta_pct_CI;
lmm_results.pval            = pval;
lmm_results.tval            = tval;
lmm_results.df              = df;
% fprintf('Group difference (A vs B): %.2f%% (95%% CI %.2f%%, %.2f%%), p = %.3g\n', ...
%         delta_pct, delta_pct_CI(1), delta_pct_CI(2), pval);
end