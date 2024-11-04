% Model y = Ax + e, E(e) = 0, E(ee') = Q

A = [1, 2;
    3, 4;
    5, 0;
    0, 6];

y = [1.5377;
    3.6948;
    -7.7193;
    7.3621];

Q = [1, 0.5, 0.5, 0.25;
    0.5, 2, 0.25, 1.0;
    0.5, 0.25, 2, 1.0;
    0.25, 1.0, 1.0, 4];

% BLUE estimate of x using 2 values of y
[x_hat_2, covariance_m_2] = BLUE_estimate(A, y, Q, 2);
x_hat_2
covariance_m_2

% BLUE estimate of x using 3 values of y
[x_hat_3, covariance_m_3] = BLUE_estimate(A, y, Q, 3);
x_hat_3
covariance_m_3

% BLUE estimate of x using all values of y
[x_hat_all, covariance_m_all] = BLUE_estimate(A, y, Q, 4);
x_hat_all
covariance_m_all

%% Compute BLUE estimate of x using n values of y
function [x_hat, covariance_m] = BLUE_estimate(A, y, Q, n)
y_n = y(1:n);
A_n = A(1:n, :);
Q_n = Q(1:n, 1:n);

x_hat   = inv(A_n' * inv(Q_n) * A_n) * A_n' * inv(Q_n) * y_n;
covariance_m = inv(A_n' * inv(Q_n) * A_n);
end