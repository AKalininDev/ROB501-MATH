%% ROB501_HW9_P1
% Computes Jacobian of a function analytically and numerically of f:R3 -> R

%% Clean Up
clear
close all
clc

%% Part A.
% Analytical Solution

% Define symbolic variables
syms x1 x2 x3
syms result

f = 3*x1*(2*x2 - x3^3) + x2^4/3;

% Compute the Jacobian of f
J = jacobian(f, [x1, x2, x3]);
disp("Analytical Jacobian:")
pretty(J)

% Evaluate the Jacobian at x_star = [1; 3; -1]
x_star = [1; 3; -1];
J_x_star = subs(J, [x1; x2; x3], x_star);
disp("Analytical Jacobian @ [1; 3; -1]:")
pretty(J_x_star)

%% Part B.
% Numerical Solution

[J_x_star_num, delta] = iterative_jacobian(@func, x_star);

disp("Numerical Jacobian @ [1; 3; -1]:")
J_x_star_num
delta

%% Part C.
% Computing the numerical Jacobian of a new function

% Verify the function
x0 = [1 2 3 4 5];
f_new = funcPartC(x0);
disp("f_new @ x0:")
f_new

% Compute the numerical Jacobian of f_new
x_star  = [1, 1, 1, 1, 1];
[J_x_star_num_new, delta_new] = iterative_jacobian(@funcPartC, x0);

disp("Numerical Jacobian @ [1, 1, 1, 1, 1]:")
J_x_star_num_new
delta_new



%% Function Definitions

% Define the function f
function f = func(x_star)
x1 = x_star(1);
x2 = x_star(2);
x3 = x_star(3);

f = 3*x1*(2*x2 - x3^3) + x2^4/3;
end


% Compute the numerical Jacobian of a function
function J_x_star_num = compute_numerical_jacobian(f, x_star, delta)

J_x_star_num = zeros(1, length(x_star));

for i = 1:length(x_star)
    x_star_plus = x_star;
    x_star_plus(i) = x_star_plus(i) + delta;
    
    x_star_minus = x_star;
    x_star_minus(i) = x_star_minus(i) - delta;
    
    f_x_star_plus = f(x_star_plus);
    f_x_star_minus = f(x_star_minus);
    
    J_x_star_num(i) = (f_x_star_plus - f_x_star_minus) / (2 * delta);
end
end

function [J_x_star, delta] = iterative_jacobian(f, x_star)
toleranceFraction = 1e-2; % 1% tolerance
delta = 1 + norm(x_star) / 10; % start with 10% of the norm of x_star

goodApprox = false;
J_x_star = compute_numerical_jacobian(f, x_star, delta);
while ~goodApprox
    delta = delta / 2;
    J_x_star_new = compute_numerical_jacobian(f, x_star, delta);
    
    difference = J_x_star_new - J_x_star;
    goodApprox = true;
    for i = 1:length(difference)
        if abs(difference(i)) > abs(J_x_star(i) * toleranceFraction)
            goodApprox = false;
        end
    end
    J_x_star = J_x_star_new;
    
end
end



