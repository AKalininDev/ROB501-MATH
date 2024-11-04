% Model y = Ax + e, E(e) = 0, E(ee') = Q

A = [1, 2;
    3, 4;
    5, 0;
    0, 6];

y = [1.5377;
    3.6948;
    -7.7193;
    7.3621];

Q_identity = eye(4);
P_100 = 100 * eye(2);
P_10E6 = 10^6 * eye(2);

% Using Deterministic Least Squares Estimation
x_hat_deterministic = LeastSquares_estimate(A, y)

% BLUE estimate of x Q_identity
[x_hat_blue, covariance_m_blue] = BLUE_estimate(A, y, Q_identity, 4)

% MVE estimate of x Q_identity and P_100
[x_hat_mve_100, covariance_m_mve_100] = MVE_estimate(A, y, Q_identity, P_100, 4)

% MVE estimate of x using Q_identity and 1E6
[x_hat_mve_1E6, covariance_m_mve_1E6] = MVE_estimate(A, y, Q_identity, P_10E6, 4)

%% Compute Least Squares Estimate Assuming Deterministic Model using Normal Equations
function x_hat = LeastSquares_estimate(A, y)
x_hat = inv(A' * A) * A' * y;
end

%% Compute BLUE estimate of x using n values of y
function [x_hat, covariance_m] = BLUE_estimate(A, y, Q, n)
y_n = y(1:n);
A_n = A(1:n, :);
Q_n = Q(1:n, 1:n);

x_hat   = inv(A_n' * inv(Q_n) * A_n) * A_n' * inv(Q_n) * y_n;
covariance_m = inv(A_n' * inv(Q_n) * A_n);
end

%% Compute MVE estimate of x using n values of y
function [x_hat, covariance_m] = MVE_estimate(A, y, Q, P, n)
y_n = y(1:n);
C_n = A(1:n, :);
Q_n = Q(1:n, 1:n);

x_hat = inv(C_n' * inv(Q_n) * C_n + inv(P)) * C_n' * inv(Q_n) * y_n;
covariance_m = P - P*C_n' * inv(C_n*P*C_n' + Q_n) * C_n * P;
end