%% ROB599_HW10_P5
% Finding R2 approximation of A

% Define the matrix A
A = [4.041,  7.046, 3.014;
    10.045, 17.032, 7.027;
    16.006, 27.005, 11.048];


% Find SVD of A
[U, S, V] = svd(A)

% Find the dA = argmin||dA||
sigma = max(S(:,end))
u_min = U(:,end)
v_min = V(:,end)
dA = -sigma * u_min * v_min'

% Find the e-value of dA'dA
e_max_val = max(eig(dA' * dA))
dA_norm = norm(dA)
A_hat = A + dA

% Verify the result
[S_hat, U_hat, V_hat] = svd(A_hat)
