%% Scatter plot example

% Example data
x = [1 2 3 4 5 6 7 8 9 10];
y = [2.1 2.8 3.5 4.2 5.4 5.9 7.1 7.8 8.6 9.4];

% Create the scatter plot
figure;
scatter(x, y, 70, 'filled');
grid on;
xlabel('X data');
ylabel('Y data');
title('Example Scatter Plot');