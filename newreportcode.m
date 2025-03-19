clear all;
close all;
clc;

% Parameters
num_users = 100;           % Number of users
num_cells = 40             % Number of small cells
area_size = 1500;         % Area size in meters (square area: 1000m x 1000m)
max_power = 1;            % Maximum transmit power per cell (Watts)
noise_power = 1e-9;       % Noise power (Watts)
bandwidth = 20e6;         % Bandwidth (20 MHz)
path_loss_exp = 3.5;      % Path loss exponent

% Randomly place users and small cells in the area
user_pos = area_size * rand(num_users, 2);    % User positions [x, y]
cell_pos = area_size * rand(num_cells, 2);    % Cell positions [x, y]

% Calculate distances between users and cells
distances = zeros(num_users, num_cells);
for u = 1:num_users
    for c = 1:num_cells
        distances(u, c) = norm(user_pos(u, :) - cell_pos(c, :));
    end
end

% Path loss model (in dB)
path_loss = 10 * path_loss_exp * log10(distances) + 20 * log10(3e8 / (4 * pi * 2.4e9));

% Received signal power (linear scale)
tx_power = max_power * ones(num_users, num_cells); % Assume full power initially
rx_power = tx_power ./ (10 .^ (path_loss / 10));

% Signal-to-Interference-plus-Noise Ratio (SINR)
sinr = zeros(num_users, num_cells);
for u = 1:num_users
    for c = 1:num_cells
        interference = sum(rx_power(u, :)) - rx_power(u, c); % Total interference
        sinr(u, c) = rx_power(u, c) / (interference + noise_power); % SINR calculation
    end
end

% Capacity (Shannon's formula)
capacity = bandwidth * log2(1 + sinr); % Capacity in bits per second

% Fairness-Based Cell Selection Mechanism
user_assignment = zeros(num_users, 1); % Store cell assignment for each user
cell_load = zeros(num_cells, 1);       % Number of users per cell

for u = 1:num_users
    % Select cell with maximum capacity, adjusted by fairness (load balancing)
    [~, sorted_cells] = sort(capacity(u, :), 'descend');
    selected_cell = -1;
    min_load = inf;
    
    % Check top 3 cells and pick the least loaded one for fairness
    for c_idx = 1:min(3, num_cells)
        cell = sorted_cells(c_idx);
        if cell_load(cell) < min_load
            min_load = cell_load(cell);
            selected_cell = cell;
        end
    end
    
    user_assignment(u) = selected_cell;
    cell_load(selected_cell) = cell_load(selected_cell) + 1;
end

% Energy Efficiency Calculation (bits per Joule)
total_power = sum(cell_load > 0) * max_power; % Total power consumed by active cells
total_capacity = sum(capacity(sub2ind(size(capacity), (1:num_users)', user_assignment)));
energy_efficiency = total_capacity / total_power;

% Fairness Index (Jain's Fairness Index)
user_rates = capacity(sub2ind(size(capacity), (1:num_users)', user_assignment));
fairness_index = (sum(user_rates)^2) / (num_users * sum(user_rates.^2));

% Results
fprintf('Total Energy Efficiency: %.2f bits/Joule\n', energy_efficiency);
fprintf('Jain''s Fairness Index: %.4f\n', fairness_index);

% User Throughput Calculation
user_throughput = user_rates; % Throughput for each user (in bits per second)

% Visualization
% Plot 1: User-Cell Assignment
figure;
scatter(cell_pos(:, 1), cell_pos(:, 2), 100, 'bs', 'filled'); % Small cells
hold on;
scatter(user_pos(:, 1), user_pos(:, 2), 50, 'ro'); % Users
for u = 1:num_users
    c = user_assignment(u);
    plot([user_pos(u, 1) cell_pos(c, 1)], [user_pos(u, 2) cell_pos(c, 2)], 'k-');
end
title('User-Cell Assignment in 5G Ultra-Dense Network');
xlabel('X (meters)');
ylabel('Y (meters)');
legend('Small Cells', 'Users'); % Only add legend for Small Cells and Users
grid on;

% Plot 2: SINR Distribution (Histogram)
figure;
histogram(sinr(:), 'Normalization', 'pdf'); % PDF of SINR values
title('SINR Distribution');
xlabel('SINR (dB)');
ylabel('Probability Density');
grid on;

% Plot 3: Capacity Distribution (Histogram)
figure;
histogram(capacity(:) / 1e6, 'Normalization', 'pdf'); % PDF of Capacity values in Mbps
title('Capacity Distribution');
xlabel('Capacity (Mbps)');
ylabel('Probability Density');
grid on;

% Plot 4: CDF of SINR (Manual Calculation)
figure;
[sorted_sinr, sinr_cdf] = calculate_cdf(sinr(:)); % Calculate CDF of SINR
plot(sorted_sinr, sinr_cdf, 'b-', 'LineWidth', 2); % Plot CDF
title('CDF of SINR');
xlabel('SINR (dB)');
ylabel('Cumulative Probability');
grid on;

% Plot 5: CDF of Capacity (Manual Calculation)
figure;
[sorted_capacity, capacity_cdf] = calculate_cdf(capacity(:) / 1e6); % Calculate CDF of Capacity
plot(sorted_capacity, capacity_cdf, 'r-', 'LineWidth', 2); % Plot CDF
title('CDF of Capacity');
xlabel('Capacity (Mbps)');
ylabel('Cumulative Probability');
grid on;

% Plot 6: User Throughput Distribution (Histogram)
figure;
histogram(user_throughput / 1e6, 'Normalization', 'pdf'); % Throughput in Mbps
title('User Throughput Distribution');
xlabel('Throughput (Mbps)');
ylabel('Probability Density');
grid on;

% Plot 7: CDF of User Throughput (Manual Calculation)
figure;
[sorted_throughput, throughput_cdf] = calculate_cdf(user_throughput / 1e6); % Calculate CDF of Throughput
plot(sorted_throughput, throughput_cdf, 'g-', 'LineWidth', 2); % Plot CDF
title('CDF of User Throughput');
xlabel('Throughput (Mbps)');
ylabel('Cumulative Probability');
grid on;

% Function to calculate CDF
function [sorted_data, cdf] = calculate_cdf(data)
    sorted_data = sort(data); % Sort the data
    cdf = (1:length(sorted_data)) / length(sorted_data); % Calculate CDF
end