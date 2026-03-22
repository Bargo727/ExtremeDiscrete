% =========================================================================
% Plot Mean Extreme First-Passage Time vs Ensemble Size (N)
% Comparison of Simulation vs. Theoretical Formula (Eq 6)
% =========================================================================

clear; clc; close all;

%% 1. SYSTEM CONFIGURATION
% -------------------------------------------------------------------------
% Using the 'Leaky Loop' model with deterministic tail (mu=1)
% to ensure finite means for all N.
d = 20;          % Minimum path length (Hard Edge)
s = 0.9;         % Probability to stay in trap (0.9 = 'sticky')

% Derived Shortest Path Probability
p_d = 1 - s;     % Probability of immediate exit

% Range of Ensemble Sizes (N) to simulate
N_list = [5, 10, 20, 50, 100, 200, 500, 1000];
num_trials = 5000; % Monte Carlo iterations per N point

%% 2. THEORETICAL CALCULATION
% -------------------------------------------------------------------------
% Formula: <T_N> = d + sum_{k=0}^{inf} exp( -N * p_d * F(k) )
% For Leaky Loop: p_d*F(k) is simply the cumulative prob P(delay <= k)
% actually p_d*F(k) = (1-s) * (1-s^(k+1))/(1-s) = 1 - s^(k+1)

% We compute the theory curve for a smooth range of N
N_smooth = logspace(log10(min(N_list)), log10(max(N_list)), 100);
theo_mean = zeros(size(N_smooth));

% Summation cutoff (when term becomes negligible)
k_max = 200; 

for i = 1:length(N_smooth)
    N_curr = N_smooth(i);
    sum_term = 0;
    
    for k = 0:k_max
        % F(k) for Leaky Loop
        F_k = (1 - s^(k+1)) / (1 - s);
        
        % The exponential approximation from your Theorem
        term = exp(-N_curr * p_d * F_k);
        sum_term = sum_term + term;
    end
    
    theo_mean(i) = d + sum_term;
end

%% 3. MONTE CARLO SIMULATION
% -------------------------------------------------------------------------
sim_mean = zeros(size(N_list));
sim_std  = zeros(size(N_list)); % For error bars

fprintf('Running Simulations...\n');
for i = 1:length(N_list)
    N = N_list(i);
    
    % Store minimum times for this N
    min_times = zeros(num_trials, 1);
    
    for t = 1:num_trials
        % 1. Simulate Delays (Geometric Dist)
        % Number of failures before first success
        delays = geornd(1-s, N, 1);
        
        % 2. Extreme Value
        T_N = d + min(delays);
        min_times(t) = T_N;
    end
    
    sim_mean(i) = mean(min_times);
    sim_std(i)  = std(min_times) / sqrt(num_trials); % Standard Error
    
    fprintf('N = %4d | <T_N> Sim: %.4f | Theory: %.4f\n', ...
            N, sim_mean(i), interp1(N_smooth, theo_mean, N));
end

%% 4. PLOTTING
% -------------------------------------------------------------------------
figure('Color','w', 'Position', [100 100 700 500]);

% 1. Plot Simulation Dots with Error Bars
errorbar(N_list, sim_mean, sim_std, 'bo', ...
    'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
hold on;

% 2. Plot Theoretical Curve
plot(N_smooth, theo_mean, 'r-', 'LineWidth', 2);

% 3. Plot the Hard Edge (d)
yline(d, 'k--', 'Hard Edge (d)', 'LabelHorizontalAlignment', 'left');

% Formatting
set(gca, 'XScale', 'log'); % Log scale for N is usually cleaner
xlabel('Ensemble Size (N)', 'FontSize', 12);
ylabel('Mean Extreme First Passage Time \langle T_N \rangle', 'FontSize', 12);
title(['Verification of \langle T_N \rangle scaling (d=' num2str(d) ', s=' num2str(s) ')']);
legend({'Simulation (Monte Carlo)', 'Theorem 1 Prediction (Eq 6)'}, 'FontSize', 11);
grid on;

% Add annotation explaining the regime
dim = [0.15 0.15 0.3 0.3];
str = {['As N \rightarrow \infty, \langle T_N \rangle \rightarrow ' num2str(d)], ...
       'Fast exponential decay predicted', 'by Entropic Penalty condition.'};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'BackgroundColor','white', 'FaceAlpha', 0.8);
