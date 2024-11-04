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
y_2 = y(1:2);
x_hat_2 = inv(A(2,2)' * inv(Q(2,2)) * A(2,2)) * A(2,2)' * inv(Q(2,2)) * y_2
covariance_m_2 = inv(A(2,2)' * inv(Q(2,2)) * A(2,2))

% BLUE estimate of x using 3 values of y
y_3 = y(1:3);
x_hat_3 = inv(A(3,2)' * inv(Q(3,3)) * A(3,2)) * A(3,2)' * inv(Q(3,3)) * y_3
covariance_m_3 = inv(A(3,2)' * inv(Q(3,3)) * A(3,2))

% BLUE estimate of x using all values of y
x_hat_all = inv(A' * inv(Q) * A) * A' * inv(Q) * y
covariance_m_all = inv(A' * inv(Q) * A)