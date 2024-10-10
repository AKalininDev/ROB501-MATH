%% ROB 501. HW#6 Problem 7 - Anatolii Kalinin
% Implementing the Matrix Inversion Lemma

%% Clear Workspace and Close All Figures

clear all; 
close all; 
clc;

%% Test the Function Computing (A + BCD)^(-1)

A = diag([0.5; 1; 1; 0.5; 1]);
inv_A = inv(A);

B = [3; 0; 2; 0; 1];
C = [0.25];
D = B';

inverse = inv_A_BCD(inv_A, B, C, D);

%% Show the Result
display(inverse)

%% Functions

%% Compute (A + BCD)^(-1)
function inverse = inv_A_BCD(inv_A, B, C, D)
    
    if ~isCompatible(inv_A, B, C, D)
        inverse = nan;

        return;
    end

    if ~isInvertible(C)
        display("C must be invertable.");
        inverse = nan;
        return;
    end

    intermediate = (inv(C) + D * inv_A*B);

    if ~isInvertible(intermediate)
        display("C^(-1) + D *A^(-1)*B) must be invertable.");
        inverse = nan;
        return;
    end

    inverse = inv_A - inv_A*B*inv(intermediate)*D*inv_A;

    return;
end

%% Check if A B C D are Compatible for the Operation
function compatible = isCompatible(inv_A, B, C, D)

    [rows_a, cols_a] = size(inv_A);
    [rows_b, cols_b] = size(B);
    [rows_c, cols_c] = size(C);
    [rows_d, cols_d] = size(D);
    

    % Check if A and C are square
    if rows_a ~= cols_a
        display("A must be a square matrix");
        compatible = false;
        return;
    end

    if rows_c ~= cols_c
        display("C must be a square matrix");
        compatible = false;
        return;
    end
    
    % Check the Dimensions of B
    if rows_b ~= rows_a
        display("num rows of B must match rows of A");
        compatible = false;
        return;
    end

    if cols_b ~= cols_c
        display("num cols of B must match cols of C");
        compatible = false;
        return;
    end

    
    % Check the Dimensions of D
    if rows_d ~= rows_c
        display("num rows of D must match rows of C");
        compatible = false;
        return;
    end

    if cols_d ~= cols_a
        display("num cols of D must match cols of A");
        compatible = false;
        return;
    end
    

    % Passed all the checks from above
    compatible = true;  
end

%% Check if Matrix is Invertable
function invertible = isInvertible(matrix)
    invertible = abs(det(matrix)) > 0;
    return;
end