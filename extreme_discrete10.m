% Comet Geometry: Extreme First-Passage Survival Probability
% Visualizing Theorem 2.1 with Grouped Bar Charts
% -------------------------------------------------------------
clear; clc; close all;

% --- Parameters ---
hub_degrees = [5 30];  % Different hub complexities
tail_length = 10;       % Distance the particle must travel unidirectionally
lambda = 2.0;           % Fixed scaling parameter lambda = N * p_d
num_sims = 100000;      % Number of single walkers to sample S(t)

max_k = 5;             % Max delay steps to visualize
k_vals = 0:max_k;

% Matrix to hold the bar chart data: rows = k values, cols = [sim_d5, sim_d30, Theory]
bar_data = zeros(length(k_vals), 3);

% --- Simulation & Theory Loop ---
for i = 1:length(hub_degrees)
    D = hub_degrees(i);
    
    % 1. Simulate single walkers (Corrected to Unidirectional Tail)
    hitting_times = zeros(num_sims, 1);
    for s = 1:num_sims
        hub_steps = 0;
        in_tail = false;
        
        while ~in_tail
            if rand() < (1/D)
                in_tail = true;
                hub_steps = hub_steps + 1; % 1 step to enter the tail
            else
                hub_steps = hub_steps + 2; % 1 step to spoke, 1 step back to hub
            end
        end
        % Total time is time spent in hub + strict unidirectional travel in tail
        hitting_times(s) = hub_steps + (tail_length - 1); 
    end
    
    % Minimum possible steps (shortest path)
    min_dist = tail_length; 
    
    % 2. Calculate single-walker arrival probability along shortest path (p_d)
    p_d = sum(hitting_times == min_dist) / num_sims;
    
    % 3. Set Ensemble Size N to keep lambda constant
    N = round(lambda / p_d);
    
    % 4. Calculate empirical P(T_N > d+k)
    for idx = 1:length(k_vals)
        k = k_vals(idx);
        S_t = sum(hitting_times > (min_dist + k)) / num_sims;
        bar_data(idx, i) = S_t^N; % Store in column 1 (d=5) or 2 (d=30)
    end
end

% --- Calculate Exact Theoretical Limit ---
% For the Comet graph, delays k only occur in even increments (2 steps per spoke visit).
% Let k = 2m. The theoretical F(k) for the Comet trap evaluates to D*(1 - (1-1/D)^(m+1)).
% We use the infinite D limit of this geometric scaling for the pure theory curve.
for idx = 1:length(k_vals)
    k = k_vals(idx);
    if mod(k, 2) ~= 0
        % Odd delays are impossible in this bipartite hub geometry, survival prob remains flat
        bar_data(idx, 3) = bar_data(idx-1, 3); 
    else
        m = k / 2;
        % Asymptotic F(k) for the injection trap
        F_k = m + 1; 
        bar_data(idx, 3) = exp(-lambda * F_k); % Store in column 3 (Theory)
    end
end

% --- Plotting the Grouped Bar Chart ---
figure('Color', 'w', 'Position', [100, 100, 800, 500]);

% Create grouped bars
b = bar(k_vals, bar_data, 'grouped', 'EdgeColor', 'k', 'LineWidth', 0.7);

% Style the bars
b(1).FaceColor = '#0072BD'; % Blue for d=5
b(2).FaceColor = '#D95319'; % Orange for d=30
b(3).FaceColor = '#EDB120'; % Yellow for Theory

% --- Formatting ---
title('Invariance of Extreme First-Passage Distributions (Comet Graph)');
xlabel('Delay steps beyond shortest path ($k$)', 'Interpreter', 'latex');
ylabel('Survival Probability $\mathbb{P}(T_N > d+k)$', 'Interpreter', 'latex');
legend({'Simulation ($d=5$)', 'Simulation ($d=30$)', 'Theoretical Limit'}, ...
       'Location', 'northeast', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);
xlim([-1, max_k+1]);
ylim([0, max(b)]);