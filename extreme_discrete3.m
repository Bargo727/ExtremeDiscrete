
clear; clc; close all;

%% 1. SYSTEM PARAMETERS
% -------------------------------------------------------------------------
% Physical parameters of the graph
s = 0.5;         % Trap 'leakiness' (Prob to stay in Head)
mu = 0.9;        % Survival probability per step in the Tail
d_val = 20;      % Total distance to target (Head + Tail)

% The Scaling Parameter (The "Knob")
lambda = 1;      % Expected number of survivors along shortest path
                 % (Keep > 1 to avoid too many empty/infinite runs)

num_trials = 50000; % Number of Monte Carlo ensembles

%% 2. CALCULATE ENSEMBLE SIZE (N)
% -------------------------------------------------------------------------
% Condition from Theorem: N * p_d -> lambda
% p_d = (1-s) * mu^(d-1)
p_d = (1 - s) * mu^(d_val - 1);

% Required N to satisfy the scaling limit
N = round(lambda / p_d);

fprintf('Running Verification for d = %d\n', d_val);
fprintf('Theoretical p_d: %.2e\n', p_d);
fprintf('Scaled Ensemble Size N: %d\n', N);
fprintf('Target Lambda: %.2f\n', lambda);
fprintf('--------------------------------------\n');

%% 3. MONTE CARLO SIMULATION
% -------------------------------------------------------------------------
min_times = zeros(num_trials, 1);
successful_trials = 0;

for i = 1:num_trials
    % 1. Simulate Delay in Head (Geometric Distribution)
    % MATLAB's geornd(p) gives number of failures (loops) before success.
    % Prob(k) = (1-p)^k * p. Here success is 'exiting', prob = (1-s).
    % So loops ~ Geometric(1-s).
    loops = geornd(1 - s, N, 1);
    
    % 2. Simulate Survival in Tail (Bernoulli Trials)
    % Each particle must survive d-1 steps with prob mu per step.
    survived_tail = rand(N, 1) < (mu^(d_val - 1));
    
    % 3. Find the Extreme Time (Minimum of Survivors)
    valid_loops = loops(survived_tail);
    
    if isempty(valid_loops)
        % No particle reached the target (Time = Infinity)
        min_times(i) = NaN; 
    else
        % T_N = distance + minimum extra loops
        min_times(i) = d_val + min(valid_loops);
        successful_trials = successful_trials + 1;
    end
end

% Filter out failed trials (where no particle arrived)
% In the scaling limit, Prob(failure) = exp(-lambda * F(inf))
arrivals = min_times(~isnan(min_times));
excess_times = arrivals - d_val; % Shift to k = 0, 1, 2...

fprintf('Simulation Complete.\n');
fprintf('Survival Rate of Ensembles: %.2f%% (Theory: %.2f%%)\n', ...
    100*length(arrivals)/num_trials, ...
    100*(1 - exp(-lambda/(1-s))) ); % Approx F(inf) for comparison

%% 4. THEORETICAL DISTRIBUTION (PMF)
% -------------------------------------------------------------------------
% F(k) for Leaky Loop: (1 - s^(k+1)) / (1 - s)
k_max = max(excess_times);
k_range = 0:k_max;

% Calculate Survival Function S(k) = P(T_N > d + k)
% Eq (5): S(k) -> exp(-lambda * F(k))
F_k = (1 - s.^(k_range + 1)) ./ (1 - s);
S_theory = exp(-lambda * F_k);

% Calculate PMF from Survival: P(T=k) = S(k-1) - S(k)
% Special case for k=0: P(T=0) = 1 - S(0)
P_theory = zeros(size(k_range));
P_theory(1) = 1 - S_theory(1); % For k=0
for j = 2:length(k_range)
    P_theory(j) = S_theory(j-1) - S_theory(j);
end

% NORMALIZE THEORY
% Since we discarded empty trials in simulation, we must normalize
% the theoretical PMF by the total probability of 'at least one arrival'.
prob_arrival_any = 1 - exp(-lambda * (1/(1-s))); % F(inf) = 1/(1-s)
P_theory_cond = P_theory / prob_arrival_any;

%% 5. MEAN VALUE VERIFICATION
% -------------------------------------------------------------------------
sim_mean = mean(arrivals);

% Theoretical Mean (Conditional)
% <T_N> = d + sum(k * P(k))
theo_mean = d_val + sum(k_range .* P_theory_cond);

fprintf('--------------------------------------\n');
fprintf('RESULTS:\n');
fprintf('Mean Extreme Time (Simulation): %.4f\n', sim_mean);
fprintf('Mean Extreme Time (Theory):     %.4f\n', theo_mean);
fprintf('Difference:                     %.4f\n', abs(sim_mean - theo_mean));

%% 6. PLOTTING
% -------------------------------------------------------------------------
figure('Color','w', 'Position', [100 100 900 400]);

% A. Histogram vs Theory
%subplot(1,2,1);
h = histogram(excess_times, 'Normalization', 'probability');
h.FaceColor = [0.2 0.4 0.8];
h.EdgeColor = 'none';
hold on;
plot(k_range, P_theory_cond, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
%title('Distribution of T_N (Excess Steps)');
xlabel('$k$','Interpreter','latex');
%ylabel('$\mathbb{P}(T_N > d + k)$','Interpreter','latex');
gg = legend('simulation', 'theory', 'Location', 'ne');
set(gg,'box','off')
set(gg,'fontsize',20)
set(gca,'fontsize',20)
set(gca, 'YScale', 'log')
axis square
grid on;
xlim([-0.5 5.5]);
ax = gca;
exportgraphics(ax, '~/Desktop/my_histogram.pdf', 'ContentType', 'Vector');

% B. Survival Function Log-Plot (The "Hard Edge" Check)
% subplot(1,2,2);
% % Empirical Survival
% [f_emp, x_emp] = ecdf(excess_times);
% plot(x_emp, 1-f_emp, 'b.', 'MarkerSize', 15); 
% hold on;
% % Theoretical Survival (Conditional)
% % S_cond(k) = S_abs(k) - P(Empty) / (1 - P(Empty)) ... roughly
% % Easier to just cumsum the conditional PMF reversed
% S_theory_cond = 1 - cumsum(P_theory_cond);
% plot(k_range, S_theory_cond, 'r-', 'LineWidth', 2);

% set(gca, 'YScale', 'log');
% %title('Survival Probability P(T_N > d+k)');
% xlabel('k (Extra Steps)');
% ylabel('log P(T_N > d+k)');
% legend('Simulation', 'Theory', 'Location', 'sw');
% grid on;
% ylim([1e-3 1]);
% xlim([0 10]);

%sgtitle(['Verification of Extreme FPT Theorem (d=' num2str(d_val) ', \lambda=' num2str(lambda) ')']);