% Verifies the Theorem in the Rare Event Limit (Poisson Regime).
% We scan Lambda = N * p_d from 0.1 to 10 to see the "transition".

clear; clc; close all;

%% 1. Setup
N = 10000;              % Fixed large N (Thermodynamic limit)
s = 0.9;                % Base 's' for geometry (defines the step size of F(k))
                        % Note: We will essentially vary the effective p_d 
                        % by treating Lambda as the control variable.

Lambda_list = logspace(-1, 1.5, 50); % Lambda from 0.1 to ~30
var_sim = zeros(size(Lambda_list));
var_theory = zeros(size(Lambda_list));

fprintf('Simulating Rare Event Regime (N=%d)...\n', N);

%% 2. Loop over Intensity (Lambda)
for i = 1:length(Lambda_list)
    Lambda = Lambda_list(i);
    
    % To simulate this Lambda with fixed N, we need a specific p_d.
    % p_d = Lambda / N.
    % In the Leaky Loop, p_d = (1-s_eff). 
    % So we tweak the "leakiness" for each trial to match the target Lambda.
    
    p_d_target = Lambda / N;
    s_eff = 1 - p_d_target; % Effective s required to get this Lambda
    
    % --- A. Monte Carlo ---
    num_trials = 50000;
    
    % Delays are Geometric(p_d_target)
    % MATLAB geornd(p) gives number of failures.
    delays = geornd(p_d_target, num_trials, N);
    min_delays = min(delays, [], 2);
    
    var_sim(i) = var(min_delays);
    
    % --- B. Theory ---
    % F(k) for geometric decay is (1 - s_eff^(k+1)) / (1 - s_eff)
    % Since s_eff is very close to 1 (because p_d is small), 
    % F(k) approximates to linear growth k+1 initially?
    % Let's use the EXACT formula from your paper Eq (13).
    
    k = 0:200;
    % F_k calculation using the s_eff for this specific point
    F_k = (1 - s_eff.^(k+1)) ./ (1 - s_eff);
    
    % Theorem 1 Distribution: P(T_N > d+k) = exp(-Lambda * F(k))
    tail_probs = exp(-Lambda * F_k);
    
    % Moments
    E_D = sum(tail_probs);
    E_D2 = sum((2.*k + 1) .* tail_probs);
    var_theory(i) = E_D2 - (E_D)^2;
end

%% 3. Plot
figure('Color','w');
semilogx(Lambda_list, var_sim, 'ko', 'MarkerFaceColor', 'b', 'DisplayName', 'Simulation');
hold on;
semilogx(Lambda_list, var_theory, 'r-', 'LineWidth', 2, 'DisplayName', 'Theory');
xlabel('Intensity \lambda = N p_d');
ylabel('Variance of T_N');
title(['Variance in the Rare Event Regime (N=' num2str(N) ')']);
legend;
grid on;

% Add annotation
text(min(Lambda_list), max(var_theory)*0.9, ...
    {'  Low \lambda: Variance is High', '  (Most particles miss target)'}, ...
    'FontSize', 10);
text(max(Lambda_list)*0.1, min(var_theory)+0.1, ...
    {'High \lambda: Variance \to 0', '(Shortest path is certain)'}, ...
    'FontSize', 10);