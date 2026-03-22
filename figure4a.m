% =========================================================================
% Entropic Collapse Visualization on the Bethe Lattice
%
% This script demonstrates that on the Bethe Lattice, the distribution of
% extreme first-passage times is NOT universal. As distance d increases,
% the survival probability drops to zero instantly, demonstrating the 
% "entropic collapse" caused by the massive weight of longer paths.
% =========================================================================
clear; clc; close all;

%% 1. PARAMETERS
% -------------------------------------------------------------------------
z = 3;                  % Coordination number
lambda_val = 1.0;       % Fixed scaling parameter (N*p_d = lambda)
num_trials = 100000;    % Monte Carlo iterations for smooth probabilities
d_values = [5,10,30]; % Distances to compare

% For a Bethe lattice, detours only occur in pairs (k = 0, 2, 4, 6...)
max_k = 5; 
k_vals = 0:2:max_k; 

% Matrix to hold the bar chart data: rows = k values, cols = [d=5, d=30, d=50]
bar_data = zeros(length(k_vals), 3);

%% 2. SIMULATION & PROBABILITY CALCULATION
% -------------------------------------------------------------------------
for i = 1:length(d_values)
    d = d_values(i);
    
    % q is the probability weight for a detour pair relative to shortest path
    q = (z - 2) / (z * (z - 1));
    
    % We consider paths up to a sufficiently large excess time
    m_vals = 0:(max_k/2 + 5); 
    
    % Expected number of particles arriving at d+2m is N * p_{d+2m}.
    % From the paper's formalism, this is approximately lambda * (q*d)^m / m!
    lambdas = lambda_val * ((q * d).^m_vals) ./ factorial(m_vals);
    
    % Generate Arrivals using Poisson distribution
    counts = poissrnd(repmat(lambdas, num_trials, 1));
    
    % Find Extreme First-Passage Time for each trial
    excess_times = zeros(num_trials, 1);
    for trial = 1:num_trials
        idx = find(counts(trial,:) > 0, 1, 'first');
        if ~isempty(idx)
            excess_times(trial) = 2 * m_vals(idx); % Convert back to k
        else
            excess_times(trial) = Inf; 
        end
    end
    
    % Calculate empirical P(T_N > d+k)
    for idx_k = 1:length(k_vals)
        k = k_vals(idx_k);
        % Survival probability is the fraction of trials where the minimum 
        % arrival time is STRICTLY GREATER than d+k
        bar_data(idx_k, i) = sum(excess_times > k) / num_trials;
    end
end

%% 3. PLOTTING THE GROUPED BAR CHART
% -------------------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on;

% Create grouped bars
b = bar(k_vals, bar_data, 'grouped', 'EdgeColor', 'k', 'LineWidth', 0.7);

% Style the bars for high contrast
b(1).FaceColor = '#0072BD'; % Blue for d=5
b(2).FaceColor = '#D95319'; % Orange for d=30
b(3).FaceColor = '#7E2F8E'; % Purple for d=50

% --- Formatting ---
title('Entropic Collapse of Extreme First-Passage Distributions (Bethe Lattice)');
xlabel('Delay steps beyond shortest path ($k$)', 'Interpreter', 'latex');
ylabel('Survival Probability $\mathbb{P}(T_N > d+k)$', 'Interpreter', 'latex');
legend({'Distance $d=5$', 'Distance $d=30$', 'Distance $d=50$'}, ...
       'Location', 'northeast', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);
xlim([-1, max_k+1]);
ylim([0, max(b)]);

% Ensure x-ticks only show the valid even values
xticks(k_vals);

hold off;
