% Model y = Ax + e, E(e) = 0, E(ee') = Q

C = [1, 2;
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

P = [0.5, 0.25;
    0.25, 0.5];

x_mean = [1; -1];
e_mean = [0; 0; 0; 0];

% BLUE estimate of x using all values of y
[x_hat_all, covariance_m_all] = MVE_estimate(A, y, Q, P, x_mean, e_mean)

%% Compute MVE estimate of x using n values of y
function [x_hat, covariance_m] = MVE_estimate(C, y, Q, P, x_mean, e_mean)
x_hat = x_mean + P*C' * inv(C*P*C' + Q) * (y - C*x_mean - e_mean);
covariance_m = P - P*C' * inv(C*P*C' + Q) * C * P;

end
