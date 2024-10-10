% ROB 501. HW#6 Problem 2A.
% Find the best fit of data and estimate its derivative

% Clear workspace and close all figures
clear all; close all; clc;

% Load data
load('DataHW06_Prob2.mat');

% Initialize variables
window_size = 3;
t_window = [];
y_window = [];
t_all = [];
y_all = [];
y_derivative = [];
y_naive_derivative = [];
last_index = 0;

% Main loop for processing data
while true
    
    [new_t, new_y, last_index] = probeNext(t, y, last_index);

    if isnan(new_t)
        break;
    end
    
    % Update windows
    t_window = insert_point(t_window, new_t, window_size);
    y_window = insert_point(y_window, new_y, window_size);
    
    % Normalize time window
    normalized_t_window = t_window - min(t_window);
    
    % Construct design matrix for cubic polynomial
    A = [ones(size(normalized_t_window')), normalized_t_window', (normalized_t_window').^2];
    
    % Find best fit coefficients
    best_fit = findBestFit(A, y_window');
    
    % Estimate derivative at the latest point
    if ~any(isnan(best_fit))
        best_fit_der = polyder(flip(best_fit));
        best_der_fit_val = polyval(best_fit_der, normalized_t_window(end));
    else
        best_der_fit_val = nan;
    end


    % Compute Naive Derivative
   
    if (size(y_window) == 1)

        old_y = y_window(1);
        old_t = t_window(1);
        y_naive_derivative(1) = 0;

    else
        new_y = y_window(end);
        new_t = t_window(end);

        y_naive_derivative(end + 1) = (new_y - old_y) / (new_t - old_t);

        old_y = new_y;
        old_t = new_t;
    
    end;
    
    
    % Store results
    y_derivative(end + 1) = best_der_fit_val;
    t_all(end + 1) = new_t;


    
end

%% Compute RMSE

valid_indices = 3:length(t);
RMSE_regression = sqrt(mean((dy(valid_indices) - y_derivative(valid_indices)).^2, "all"));
RMSE_naive = sqrt(mean((dy(valid_indices) - y_naive_derivative(valid_indices)).^2, "all"));


%% Plot the results
figure('Position', [100, 100, 1600, 900]);
sgtitle('Figure 1. Problem 2, Comparing dy/dt Estimation Techniques',...
    'FontSize', 24, 'FontWeight', 'bold');

% Show Original Data
subplot(3,1,1);

original_data_name = "Original Data.";
plot(t, y, 'b--', 'LineWidth', 4, "DisplayName", original_data_name);
grid on;
set(gca, 'FontSize', 14);
set(gca,'Xticklabel',[]);
ylabel('y', 'FontSize', 16);
legend(FontSize=16);

% Naive Estimation vs True Derivative
subplot(3,1,2);

% naive data
naive_der_est_name = "Naive dy/dt Estimation. (Part A)";
plot(t_all(3:end), y_naive_derivative(3:end), 'b-', 'LineWidth', 4, 'DisplayName', naive_der_est_name);
hold on;

% true derivative
true_der_name = "True Derivative.";
plot(t, dy, 'r--', 'LineWidth', 3, "DisplayName", true_der_name);
hold off;


% decorations
grid on;
set(gca, 'FontSize', 14);
ylabel('dy/dt', 'FontSize', 16);
set(gca,'Xticklabel',[]);
legend(FontSize=16);
x_pos = 0.01;  
y_pos = 10; 
text(x_pos, y_pos, sprintf('$\\sqrt{\\frac{1}{n}\\sum_{i=1}^n (\\hat{y}_i - y_i)^2}$ = %.4f', RMSE_naive), ...
    'Interpreter', 'latex', 'FontSize', 14, 'BackgroundColor', 'white', 'EdgeColor', 'black');



% Regression Derivative vs Real Derivative
subplot(3,1,3);

% Regression Estimation
estimated_der_name = "Regression dy/dt Estimation. (Part B)";
plot(t_all(3:end), y_derivative(3:end), 'b-', 'LineWidth', 4, 'DisplayName', estimated_der_name);
hold on;

true_der_name = "True Derivative.";
plot(t, dy, 'r--', 'LineWidth', 3, "DisplayName", true_der_name);
hold off;
grid on;
set(gca, 'FontSize', 14);
xlabel('t', 'FontSize', 16);
ylabel('dy/dt', 'FontSize', 16);
legend(FontSize=16);
x_pos = 0.01;  
y_pos = 10; 
text(x_pos, y_pos, sprintf('$\\sqrt{\\frac{1}{n}\\sum_{i=1}^n (\\hat{y}_i - y_i)^2}$ = %.4f', RMSE_regression), ...
    'Interpreter', 'latex', 'FontSize', 14, 'BackgroundColor', 'white', 'EdgeColor', 'black');



% Adjust layout
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

% Save the figure
print('ROB501-HW6-P2-Derivatives.png', '-dpng', '-r300');

% Function to add new data and replace old data in the moving window
function window = insert_point(window, new_point, window_size)
    if length(window) < window_size
        window(end + 1) = new_point;
    else
        window = [window(2:end), new_point];
    end
end

% Function to find best fit coefficients
function best_fit_coeff = findBestFit(A, y)
    
    if (det(A'*A) ~= 0)
        best_fit_coeff = inv(A'*A)*A'* y;
    else
        best_fit_coeff = nan(size(A,2));
    end


end

% Function to get next data point
function [new_t, new_y, last_index] = probeNext(t, y, last_index)
    
    last_index = last_index + 1;


    if last_index <= length(t)
        new_t = t(last_index);
        new_y = y(last_index);
    else
        new_t = nan;
        new_y = nan;
    end
end