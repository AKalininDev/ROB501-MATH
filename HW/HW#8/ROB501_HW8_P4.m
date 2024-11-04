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

% MVE estimate of x using 1 value of y
[x_hat_1, covariance_m_1] = MVE_estimate(A, y, Q, P, 1)

% MVE estimate of x using 2 values of y
[x_hat_2, covariance_m_2] = MVE_estimate(A, y, Q, P, 2)

% BLUE estimate of x using 3 values of y
[x_hat_3, covariance_m_3] = MVE_estimate(A, y, Q, P, 3)

% BLUE estimate of x using all values of y
[x_hat_all, covariance_m_all] = MVE_estimate(A, y, Q, P, 4)

%% Compute MVE estimate of x using n values of y
function [x_hat, covariance_m] = MVE_estimate(A, y, Q, P, n)
y_n = y(1:n);
C_n = A(1:n, :);
Q_n = Q(1:n, 1:n);

x_hat = inv(C_n' * inv(Q_n) * C_n + inv(P)) * C_n' * inv(Q_n) * y_n;
covariance_m = P - P*C_n' * inv(C_n*P*C_n' + Q_n) * C_n * P;

end
