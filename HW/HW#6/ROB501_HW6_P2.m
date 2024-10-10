%% ROB 501. HW#6 Problem 2.
% Test Naive and Regression Derivative estimation Techniques

%% Clear Workspace and Close All Figures

clear all; 
close all; 
clc;

%% Load data & Set the Simulation Parameters

load('DataHW06_Prob2.mat'); % Loads t, y, dy
window_size = 3;

%% Process Incoming Data

[t_processed, y_derivative, y_naive] = processData(t, y, window_size);

RMSE_regression = computeRMSE(dy(window_size:end), y_derivative'); % Need to make sure that dy and the y_derivative are the same size
RMSE_naive = computeRMSE(dy(window_size:end), y_naive');

%% Plot the Results

figure('Position', [100, 100, 1600, 900]);
sgtitle('Figure 1. Problem 2, Comparing dy/dt Estimation Techniques',...
    'FontSize', 24, 'FontWeight', 'bold');

% Show Original Data
subplot(3,1,1);

original_data_name = "Original Data.";
plot(t, y, '*', 'MarkerSize', 4, "DisplayName", original_data_name);

grid on;
set(gca, 'FontSize', 14);
set(gca,'Xticklabel',[]);
ylabel('y', 'FontSize', 16);
legend(FontSize=16);

% Naive Estimation vs True Derivative
subplot(3,1,2);

naive_der_est_name = "Naive dy/dt Estimation. (Part A)";
plot(t_processed, y_naive, 'b-', 'LineWidth', 4, 'DisplayName', naive_der_est_name);
hold on;

true_der_name = "True Derivative.";
plot(t, dy, 'r--', 'LineWidth', 3, "DisplayName", true_der_name);
hold off;

x_pos = 0.01;  
y_pos = 10; 
text(x_pos, y_pos, sprintf('$\\sqrt{\\frac{1}{n}\\sum_{i=1}^n (\\hat{y}_i - y_i)^2}$ = %.4f', RMSE_naive), ...
    'Interpreter', 'latex', 'FontSize', 14, 'BackgroundColor', 'white', 'EdgeColor', 'black');

grid on;
set(gca, 'FontSize', 14);
ylabel('dy/dt', 'FontSize', 16);
set(gca,'Xticklabel',[]);
legend(FontSize=16);


% Regression Derivative vs Real Derivative
subplot(3,1,3);

estimated_der_name = "Regression (3 points) dy/dt Estimation. (Part B)";
plot(t_processed, y_derivative, 'b-', 'LineWidth', 4, 'DisplayName', estimated_der_name);
hold on;

true_der_name = "True Derivative.";
plot(t, dy, 'r--', 'LineWidth', 3, "DisplayName", true_der_name);
hold off;

x_pos = 0.01;  
y_pos = 10; 
text(x_pos, y_pos, sprintf('$\\sqrt{\\frac{1}{n}\\sum_{i=1}^n (\\hat{y}_i - y_i)^2}$ = %.4f', RMSE_regression), ...
    'Interpreter', 'latex', 'FontSize', 14, 'BackgroundColor', 'white', 'EdgeColor', 'black');

grid on;
set(gca, 'FontSize', 14);
xlabel('t', 'FontSize', 16);
ylabel('dy/dt', 'FontSize', 16);
legend(FontSize=16);

% Save the figure
print('ROB501-HW6-P2-Derivatives.png', '-dpng', '-r300');

%% Various Functions

%% Function that Simulates Data Postprocessing in Real Time
function [t_processed, y_derivative, y_naive] = processData(t, y, window_size)
     
    % Output Arrays
    y_derivative = [];
    y_naive = [];
    t_processed = [];

    % Define Arrays for Moving Window
    t_window = [];
    y_window = [];

    % Temp Variables
    last_index = 0;
    old_y = NaN;
    old_t = NaN;
    new_y = NaN;
    new_t = NaN;

    while (~isnan(new_y)) || (last_index == 0)

        % Update old and new values
        old_t = new_t;
        old_y = new_y;
        [new_t, new_y, last_index] = probeNext(t, y, last_index);

        if isnan(new_y)
            break
        end

        % Update windows
        t_window = insert_point(t_window, new_t, window_size);
        y_window = insert_point(y_window, new_y, window_size);
            
        % Compute Naive and Regression Derivative
        naive_der_val = estimateNaiveDerivative(new_t, new_y, old_t, old_y);
        
        % Compute Second Order Regression Derivative
        regression_der_val = estimateSecondOrderRegressionDerivative(t_window, y_window);

        % Store results
        y_naive(end + 1) = naive_der_val;
        y_derivative(end + 1) = regression_der_val;
    
    end

    % Remove the initial results, where there is not enough data to estimate the derivatve.
    t_processed = t(window_size:end);
    y_naive = y_naive(window_size:end);
    y_derivative = y_derivative(window_size:end);

end

%% Simulate Probing Next Datapoint
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

%% Insert New Point into a Moving Window of Fixed Size
function window = insert_point(window, new_point, window_size)
    if length(window) < window_size
        window(end + 1) = new_point;
    else
        window = [window(2:end), new_point];
    end
end

%% Estimate Naive Derivative
function naive_der_val = estimateNaiveDerivative(new_t, new_y, old_t, old_y)
    if isnan(old_t) || isnan(old_y)
        naive_der_val = 0;
    else
        naive_der_val = (new_y - old_y) / (new_t - old_t);
    end
end

%% Estimate Second Order Regression Derivative
function best_der_fit_val = estimateSecondOrderRegressionDerivative(t_window, y_window)
    
    normalized_t_window = t_window - min(t_window);
      
    A = [ones(size(normalized_t_window')), normalized_t_window', (normalized_t_window').^2];
    best_fit = findBestFit(A, y_window');

    if ~any(isnan(best_fit))
        best_fit_der = polyder(flip(best_fit));
        best_der_fit_val = polyval(best_fit_der, normalized_t_window(end));
    else
        best_der_fit_val = nan;
    end

end

%% Find Best Fit Regression
function best_fit_coeff = findBestFit(A, y)
    
    if (det(A'*A) ~= 0)
        best_fit_coeff = inv(A'*A)*A'* y;
    else
        best_fit_coeff = nan(size(A,2));
    end
end

%% Compute RMSE
function RMSE = computeRMSE(trueVals, estimatedVals)
    RMSE = sqrt(mean((trueVals - estimatedVals).^2, "all"));
    
end
