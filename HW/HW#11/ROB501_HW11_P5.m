%% ROB501_HW11_P5
% Implementing Newton-Raphson Method to solve F(x) = 0

%% Define the Initial Guess and Solver Parameters
x_initial = [10; 10];
eps = 1e-10; % Tolerance for convergence


%% Solve using Newton-Raphson Method
x_current = x_initial;
d_norm = inf;

% Newton-Raphson iteration
while d_norm > eps
    try
        delta = dF(x_current) \ F(x_current);
    catch
        disp('Jacobian is singular, cannot proceed with Newton-Raphson method.');
        break;
    end
    d_norm = norm(delta);
    x_current = x_current - delta;
end

%% Show the Results
x_initial
x_current
F(x_current)


%% Function Definitions

% Function to compute F(x)
function F_val = F(x)
x1 = x(1);
x2 = x(2);
F_val = [3 + 4*x1 + 3*x2 - x1^2 - 2*x1*x2;
    4 + 2*x1 + x2 - x1*x2 - 2*x2^2];
end

% Function to compute the Jacobian matrix of F(x)
function dF_val = dF(x)
x1 = x(1);
x2 = x(2);
dF_val = [4 - 2*x1 - 2*x2, 3 - 2*x1;
    2 - x2, 1 - x1 - 4*x2];
end
