%% ROB501-HW#8-P6
% Computes integrals and gram matricies for minimum norm solution.

%% Clean Up
clear
close all
clc

%% Define Matrix A for QR Factorization

A = [1, 2;
    3, 4;
    5, 6];

%% Computing QR Factorization "by hand"

v1 = A(:,1);
q1 = v1/norm(v1);

v2 = A(:,2) - (A(:,2)' * q1) * q1;
q2 = v2/norm(v2);

Q = [q1, q2]

r11 = A(:, 1)' * q1;
r21 = A(:, 1)' * q2;
r12 = A(:, 2)' * q1;
r22 = A(:, 2)' * q2;

R = [r11, r12;
    r21, r22]


%% Using MATLAB's QR Factorization

[Q_Mreg, R_Mreg] = qr(A)

%% Using MATLAB's econ QR Factorization

[Q_Mecon, R_Mecon] = qr(A, 0)