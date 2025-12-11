function [idx_sample, y, y_c,y_o, probabilities] = sampling_vonMises(theta, b_PF, tao_cardinal, tao_oblique, nSample)
if b_PF >= 0
    %%% overrepresent cardinal orientation
    y =  (exp(b_PF * cos(2 * theta))+exp(b_PF * cos(2*(theta - pi/2))));
else
    %%%% overrepresent oblique orientation
    y =  (exp(b_PF * cos(2 * (theta - pi/4)))+exp(b_PF * cos(2*(theta - 3 *pi/4))));
end

if tao_cardinal >= 0
    %%%% overrepresent "positive" cardinal orientation (0 degree)
    y_c =  (exp(tao_cardinal * cos(theta)) + exp(tao_cardinal * cos( theta - pi)));
else
    %%%% overrepresent "negative" cardinal orientation (90 degree)
    y_c =  (exp(tao_cardinal * cos(theta + pi/2)) + exp(tao_cardinal * cos( theta - pi/2)));
end

if tao_oblique >= 0
    %%%% overrepresent "positive" oblique orientation (45 degree)
    y_o =  (exp(tao_oblique * cos(theta - pi/4)) + exp(tao_oblique * cos( theta + 3 * pi/4)));
else
    %%%% overrepresent "negative" oblique orientation (135 degree)
    y_o =  (exp(tao_oblique * cos(theta - 3 * pi/4)) + exp(tao_oblique * cos( theta + pi/4)));
end

y = y / sum(y);
y_c = y_c / sum(y_c);
y_o = y_o / sum(y_o);

y_all = y .* y_c .* y_o;

y_all = y_all / sum(y_all); % Normalize to sum to 1
probabilities = y_all; 

% Create cumulative distribution function (CDF)
cumdf =  cumsum(probabilities);
cumdf(1) = 0; % Add 0 at the start

% Generate uniform random samples
linearProb = linspace(0,0.99999999,nSample);
phi        = interp1(cumdf,theta, linearProb, 'nearest');

idx_sample = arrayfun(@(x)find(theta == x), phi);
idx_sample = unique(idx_sample);
    

end
