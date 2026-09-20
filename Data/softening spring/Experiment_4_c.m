%% MECH 6140 - HW1
% Section 4c - Softening spring representative cases

clear;
clc;
close all;

%% Softening spring
fint = @(u) 12*u - u.^3;

uCurve = linspace(-5,5,1000);
PCurve = fint(uCurve);

files = [
    "Softening_spring_P=12_1_load_step.mat"
    "Softening_spring_P=20_1_load_step.mat"
];

labels = [
    "Convergent case: P = 12, u_0 = 0"
    "Unintended equilibrium: P = 20, u_0 = 0"
];

colors = [
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
];

for i = 1:2

    data = load(files(i), 'results', 'options');

    r = data.results.exact;
    P = data.options.finalLoad;

    uIter = r.history.displacement;
    pIter = r.history.targetLoad;

    valid = isfinite(uIter) & isfinite(pIter);

    figure;
    hold on;

    plot(uCurve, PCurve, 'k-', 'LineWidth', 2, ...
        'DisplayName', 'Softening spring');

    yline(P, '--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('P = %.0f', P));

    plot(uIter(valid), pIter(valid), '-o', ...
        'Color', colors(i,:), ...
        'LineWidth', 1.8, ...
        'MarkerSize', 6, ...
        'DisplayName', 'Newton iterations');

    xlabel('Displacement, u', 'FontWeight','bold');
    ylabel('Load, P', 'FontWeight','bold');

    title(labels(i), 'FontWeight','normal');

    ylim([-30 30]);

    legend('Location','best');
    grid on;
    box on;
    hold off;

end