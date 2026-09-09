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

%% Scatter plot with a different color for each series

% Each column contains one scatter series
ySeries = [
	2.1  2.4  1.8;
	2.8  3.1  2.5;
	3.5  3.8  3.2;
	4.2  4.5  3.9;
	5.4  5.1  4.8;
	5.9  6.4  5.5;
	7.1  6.8  6.5;
	7.8  8.2  7.3;
	8.6  8.9  8.1;
	9.4  9.7  8.8
];

numberOfSeries = size(ySeries, 2);
seriesColors = lines(numberOfSeries);

figure;
hold on;
for seriesIndex = 1:numberOfSeries
	scatter(x, ySeries(:, seriesIndex), 70, seriesColors(seriesIndex, :), 'filled');
end
hold off;
grid on;
xlabel('X data');
ylabel('Y data');
title('Scatter Plot with Different Colors');
legend('Series 1', 'Series 2', 'Series 3', 'Location', 'best');

%% Weird scatter pattern

numberOfPoints = 250;
theta = linspace(0, 8 * pi, numberOfPoints);
radius = 1.5 + 0.8 * sin(5 * theta) + 0.25 * cos(13 * theta);

weirdX = radius .* cos(theta);
weirdY = radius .* sin(theta);
markerSizes = 20 + 80 * (0.5 + 0.5 * sin(3 * theta));

figure;
scatter(weirdX, weirdY, markerSizes, theta, 'filled');
axis equal;
grid on;
colormap turbo;
colorbar;
xlabel('X position');
ylabel('Y position');
title('Weird Warped Star Scatter Pattern');