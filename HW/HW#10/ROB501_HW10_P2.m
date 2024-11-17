%% ROB501_HW10_P2
% Implement Kalman Filter using provided data

%% Clean Up
clear
close all
clc

%% Part A.

% Load data
load KFHW10data/SegwayData4KF.mat

% x = [phi, theta, phi_dot, theta_dot]

% Run Kalman Filter
x_k = x0;
P_k = P0;
K = zeros(4, length(t));

for i = 1:length(t)
    % Update Kalman Filter
    [K_k, x_k, P_k] = updateKF(A, B, C, x_k, P_k, Q, G, R, u(i), y(i));
    
    % Store data
    x(:,i) = x_k;
    
    % Store Kalman Gain
    K(:,i) = K_k;
    
end

%% Part B
% Plot Results

% Define colors
phi_color = [0.3010, 0.7450, 0.9330];
theta_color = [0.8500, 0.3250, 0.0980];
k1_color = [0.4660, 0.6740, 0.1880];
k2_color = [0.8500, 0.3250, 0.0980];
k3_color = [0.9290, 0.6940, 0.1250];
k4_color = [0.4940, 0.1840, 0.5560];

% Plot phi, theta vs time on first plot
figure('Position', [0, 0, 1200, 1000]);
tLayout = tiledlayout(3,1,'Padding','Compact');
tLayout.Title.String = "Kalman Filter Results";
tLayout.Title.FontSize = 20;
tLayout.Title.FontWeight = 'bold';

% Plot phi and theta
nexttile;
plot(t, x(1,:), 'Color', phi_color, 'LineWidth', 3, 'DisplayName', '$\phi$');
hold on;
plot(t, x(2,:), 'Color', theta_color, 'LineWidth', 3, 'DisplayName', '$\theta$');
hold off;
setSubplotProperties(gca);
ylabel('$\phi$, $\theta$ (rad)', 'FontSize', 18, 'Interpreter', 'latex');

% Plot phi_dot, theta_dot
nexttile;
plot(t, x(3,:), 'Color', phi_color, 'LineWidth', 3, 'DisplayName', '$\dot{\phi}$');
hold on;
plot(t, x(4,:), 'Color', theta_color, 'LineWidth', 3, 'DisplayName', '$\dot{\theta}$');
hold off;
setSubplotProperties(gca);
ylabel('$\dot{\phi}$, $\dot{\theta}$ (rad/s)', 'FontSize', 18, 'Interpreter', 'latex');

% Plot the four components of the Kalman Gain
nexttile;
plot(t, K(1,:), 'Color', k1_color, 'LineWidth', 3, 'DisplayName', '$K_1$');
hold on;
plot(t, K(2,:), 'Color', k2_color, 'LineWidth', 3, 'DisplayName', '$K_2$');
plot(t, K(3,:), 'Color', k3_color, 'LineWidth', 3, 'DisplayName', '$K_3$');
plot(t, K(4,:), 'Color', k4_color, 'LineWidth', 3, 'DisplayName', '$K_4$');
hold off;
setSubplotProperties(gca);
ylabel('Kalman Gains', 'FontSize', 18, 'Interpreter', 'latex');

print('ROB501-HW10-P2-Results.png', '-dpng', '-r300');

%% Part C
% Report the values to which the Kalman gains converge. Compare to Kss and Pss from dlqe()

% Compute the final average values of the Kalman Gain
K_last = K(:, end)

% Compute steady state Kalman Gain and Covariance
[Kss, Pss] = dlqe(A, G, C, R, Q)


%% Helper Function
function setSubplotProperties(ax)
grid(ax, 'on');
grid(ax, 'minor');
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 16;
xlabel(ax, '\textbf{Time (s)}', 'FontSize', 20, 'Interpreter', 'latex');
legend(ax, 'FontSize', 18, 'Location', 'northeast', 'Interpreter', 'latex');
end



function [K_k, x_k, P_k] = updateKF(A, B, C, x, P, Q, G, R, u, y)
% Update Kalman Filter
% Inputs:
%   A: State matrix
%   B: Input matrix
%   C: Output matrix
%   D: Feedthrough matrix
%   G: Disturbance matrix
%   x: State estimate
%   P: Covariance estimate
%   Q: Process noise covariance
%   R: Measurement noise covariance
%   u: Input
%   y: Measurement
% Outputs:
%   K: Kalman gain
%   x: Updated state estimate
%   P: Updated covariance estimate

% Prediction

K_k = P*C'*inv(C*P*C' + Q);
x_k = A*x + B*u;
x_k = x_k + A*K_k*(y - C*x);

P_k = A*(P - K_k*C*P)*A' + G*R*G';

end