%% ROB501-HW#8-P3
% Computes integrals and gram matricies for minimum norm solution.

%% Clean Up
clear
close all
clc

%% Define the Basis for X
syms t;

y1 = 1;
y2 = t;
y3 = t^2;
y4 = sin(pi*t);

%% Part A

% Solving <f, y2> = 2, f = argmin ||f||

Gt = inner_product(y2, y2)
b = 2 / Gt;
f = b * y2

%% Part B

% Solving <f, y2> = 1, f = argmin ||f||
% Solving <f, y4> = pi, f = argmin ||f||

Gt = [inner_product(y2, y2), inner_product(y4, y2);
    inner_product(y2, y4), inner_product(y4, y4)]
b = inv(Gt)*[2; pi];

f = vpa(b(1) * y2 + b(2) * y4, 3)

%% Define Inner Product

function [ip] = inner_product(f, g)
syms t;
ip = int(f*g, t, 0, 2);
ip = eval(ip);
end
