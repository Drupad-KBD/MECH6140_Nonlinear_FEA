%% MECH 6140 - HW1
% Sections 4a and 4b

clear;
clc;
close all;

%% ================================================================
% COMMON PLOT SETTINGS
% ================================================================

methods = ["exact", "modified", "secant", "numerical"];

colors = [
    0.0000  0.4470  0.7410;   % exact     - blue
    0.8500  0.3250  0.0980;   % modified  - red
    0.9290  0.6940  0.1250;   % secant    - orange
    0.4660  0.6740  0.1880];  % numerical - green

% Same line-width pattern for all plots
lineWidths = [3.0, 2.0, 2.0, 1.5];


%% ================================================================
% SECTION 4a - TANGENT METHOD COMPARISON
% ================================================================

data = load( ...
    'HW1_NonlinearEquilibrium_Results_1_load_step.mat', ...
    'results');

results = data.results;


%% Plot 1 - Full residual range

figure;
hold on;

for i = 1:numel(methods)

    name = methods(i);

    iteration = results.(name).history.iteration;
    Ri = abs(results.(name).history.residual);

    valid = isfinite(Ri) & Ri > 0;

    plot(iteration(valid), Ri(valid), ...
        '-', ...
        'Color', colors(i,:), ...
        'LineWidth', lineWidths(i), ...
        'DisplayName', char(name));
end

xlabel('Iteration number','FontWeight','bold');
ylabel('|R_i|','FontWeight','bold');

title('Residual Comparison','FontWeight','normal');

legend('Location','best');

grid on;
box on;
hold off;


%% Plot 2 - Zoomed residual range

figure;
hold on;

for i = 1:numel(methods)

    name = methods(i);

    iteration = results.(name).history.iteration;
    Ri = abs(results.(name).history.residual);

    valid = isfinite(Ri) & Ri > 0;

    plot(iteration(valid), Ri(valid), ...
        '-', ...
        'Color', colors(i,:), ...
        'LineWidth', lineWidths(i), ...
        'DisplayName', char(name));
end


% Zoom based only on converged methods
zoomMethods = ["exact", "secant", "numerical"];

allIterations = [];
allResiduals = [];

for i = 1:numel(zoomMethods)

    name = zoomMethods(i);

    iteration = results.(name).history.iteration;
    Ri = abs(results.(name).history.residual);

    valid = isfinite(Ri) & Ri > 0;

    allIterations = [allIterations; iteration(valid)];
    allResiduals = [allResiduals; Ri(valid)];
end

xlim([min(allIterations) max(allIterations)]);
ylim([0 1.15*max(allResiduals)]);

xlabel('Iteration number','FontWeight','bold');
ylabel('|R_i|','FontWeight','bold');

title('Residual Comparison - Zoomed','FontWeight','normal');

legend('Location','best');

grid on;
box on;
hold off;


%% Tangent method results

fprintf('\nTangent Method Results\n');
fprintf('-----------------------------\n');

for i = 1:numel(methods)

    name = methods(i);

    if results.(name).converged

        fprintf('%-10s : %d iterations\n', ...
            char(name), ...
            results.(name).numberOfCorrections);

    else

        fprintf('%-10s : Failed to converge\n', ...
            char(name));
    end
end


%% ================================================================
% SECTION 4b - LOAD STEPPING
% ================================================================

loadSteps = [1 2 5 10];

files = [
    "HW1_NonlinearEquilibrium_Results_1_load_step.mat"
    "HW1_NonlinearEquilibrium_Results_2_load_step.mat"
    "HW1_NonlinearEquilibrium_Results_5_load_step.mat"
    "HW1_NonlinearEquilibrium_Results_10_load_step.mat"
];

totalIterations = NaN(numel(methods), numel(loadSteps));

for i = 1:numel(loadSteps)

    data = load(files(i), 'results');

    for j = 1:numel(methods)

        name = methods(j);

        if data.results.(name).converged

            totalIterations(j,i) = ...
                data.results.(name).numberOfCorrections;
        end
    end
end


%% Plot 3 - Exact method only

figure;

plot(loadSteps, totalIterations(1,:), ...
    '-', ...
    'Color', colors(1,:), ...
    'LineWidth', lineWidths(1));

xlabel('Number of load increments','FontWeight','bold');
ylabel('Total Newton iterations','FontWeight','bold');

title('Load Stepping','FontWeight','normal');

xticks(loadSteps);

grid on;
box on;


%% Plot 4 - All methods

figure;
hold on;

for i = 1:numel(methods)

    plot(loadSteps, totalIterations(i,:), ...
        '-', ...
        'Color', colors(i,:), ...
        'LineWidth', lineWidths(i), ...
        'DisplayName', char(methods(i)));
end

xlabel('Number of load increments','FontWeight','bold');
ylabel('Total Newton iterations','FontWeight','bold');

title('Load Stepping Comparison','FontWeight','normal');

xticks(loadSteps);

legend('Location','best');

grid on;
box on;
hold off;


%% Load stepping results

fprintf('\nLoad Stepping Results\n');
fprintf('-------------------------------------------------\n');
fprintf('Steps       Exact     Modified     Secant     Numerical\n');

for i = 1:numel(loadSteps)

    fprintf('%-11d %-9g %-12g %-10g %-10g\n', ...
        loadSteps(i), ...
        totalIterations(1,i), ...
        totalIterations(2,i), ...
        totalIterations(3,i), ...
        totalIterations(4,i));
end