%% ROB 501. HW#5 Problem 4C.
% Find the best fit (sinusoidal functions) of data

%% Load the Data
load HW05_Prob4_Data.mat

%% Define Y, A
y = f;

A = [sin(pi *t), sin(2*pi *t), sin(3*pi *t), sin(4*pi *t), sin(5*pi *t)]

%% Check for 0 determinant

if (det(A'*A) ~= 0)
    display("det(A^TA) != 0. Can compute the best fit.")
end

%% Compute the Coefficients
a = inv(A'*A)*A'* y;

% Convert to standard MATLAB representation
p = flipud(a)

%% Plot the Results on the same Grid

% Evaluate the function fit in the range of interest
t_values = linspace(min(t) - 1.5, max(t) + 1.5, 100000); 
y_values = a(1)*sin(pi *t_values) + a(2)*sin(2*pi *t_values) + a(3)*sin(3*pi *t_values) + ...
           a(4)*sin(4*pi *t_values) + a(5)*sin(5*pi *t_values);

% Generate the figure
figure('Position', [100, 100, 1600, 900]); 

% Plot the original data and the fitted line
hold on;
plot(t, y, 's', "MarkerSize", 15, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');  % Original data (blue markers)
plot(t_values, y_values, "LineWidth", 5);
hold off;

% Add grid
grid("on");

% Adjust tick size
set(gca, 'FontSize', 18);

% Add axis labels
xlabel('t', 'FontSize', 20);
ylabel('y', 'FontSize', 20);

% Generate the equation as a string for the legend
coeff_str = sprintf('%.2f \\sin(5\\pi t)', p(1));  % Start with the highest order term

for i = 2:length(p)
    power = 5 - (i - 1);
    if p(i) >= 0
        coeff_str = sprintf('%s + %.2f \\sin(%d\\pi * t)', coeff_str, p(i), power);
    else
        coeff_str = sprintf('%s - %.2f \\sin(%d\\pi * t)', coeff_str, abs(p(i)), power);
    end
end

% Display the equation in the legend using LaTeX
legend_str = sprintf('Best Fit (MSE): p(t)= $%s$', coeff_str);
legend({'Original Data', legend_str}, 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

title('Figure 2. Problem 4C. Least Squares Sinusoidal Fit.', 'FontSize', 28, 'FontWeight', 'bold');

% Reduce Margins
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02));

% Save the figure
print('ROB501-HW#4-P4C.png', '-dpng', '-r300');