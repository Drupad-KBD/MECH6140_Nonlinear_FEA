%% MECH 6140 - HW1
% Section 3 - VERIFY

clear;
clc;
close all;

%% Problem definition

fint = @(u) 10*u + 10*u.^3;
KtExact = @(u) 10 + 30*u.^2;

P = 20;
uKnown = 1;

%% Load saved exact-Newton result

load('HW1_NonlinearEquilibrium_Results.mat', 'results', 'options');

exactResult = results.exact;


%% ================================================================
% VERIFY 1 - Known equilibrium solution
% ================================================================

uNumerical = exactResult.finalDisplacement;
uError = abs(uNumerical - uKnown);

fprintf('\nVERIFY 1\n');
fprintf('Known u     = %.10f\n', uKnown);
fprintf('Numerical u = %.10f\n', uNumerical);
fprintf('Error       = %.3e\n', uError);

uPlot = linspace(0, 1.5, 400);

figure;
plot(uPlot, fint(uPlot), 'LineWidth', 1.5);
hold on;
yline(P, '--', 'LineWidth', 1.2);
plot(uKnown, P, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(uNumerical, P, 'x', 'MarkerSize', 9, 'LineWidth', 1.5);

formatPlot( ...
    'Displacement, u', ...
    'Force', ...
    'Verification of Hardening-Spring Equilibrium');

legend('Internal force', ...
       'Applied load', ...
       'Known solution', ...
       'Numerical solution', ...
       'Location','best');

hold off;


%% ================================================================
% VERIFY 2 - Analytical vs numerical tangent
% ================================================================

uCheck = 1;

KtAnalytical = KtExact(uCheck);

hValues = logspace(-16, -1, 200);

KtNumerical = zeros(size(hValues));
relativeError = zeros(size(hValues));

for i = 1:length(hValues)

    KtNumerical(i) = computeNumericalTangent( ...
        uCheck, fint, hValues(i));

    relativeError(i) = ...
        abs(KtNumerical(i) - KtAnalytical) / abs(KtAnalytical);
end

valid = isfinite(relativeError) & relativeError > 0;

fprintf('\nVERIFY 2\n');
fprintf('Analytical tangent at u = 1 = %.10f\n', KtAnalytical);


% Analytical vs numerical tangent

figure;
semilogx(hValues, KtNumerical, 'LineWidth', 1.5);
hold on;
yline(KtAnalytical, '--', ...
    'Analytical tangent', ...
    'LineWidth', 1.2);

formatPlot( ...
    'Scalar perturbation size, h', ...
    'Tangent stiffness', ...
    'Analytical and Numerical Tangent Comparison at u = 1');

legend('Numerical tangent', ...
       'Analytical tangent', ...
       'Location','best');

hold off;


% Relative tangent error

figure;
loglog( ...
    hValues(valid), ...
    relativeError(valid), ...
    'o-', ...
    'MarkerSize', 4, ...
    'LineWidth', 1.2);

formatPlot( ...
    'Scalar perturbation size, h', ...
    'Relative tangent error', ...
    'Relative Numerical Tangent Error at u = 1');


%% ================================================================
% VERIFY 3 - Ri versus iteration number
% ================================================================

iteration = exactResult.history.iteration;
Ri = abs(exactResult.history.residual);

fprintf('\nVERIFY 3\n');
fprintf('Residual tolerance = %.3e\n', options.residualTolerance);

figure;
plot(iteration, Ri, 'o-', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.5);

hold on;
yline(0, '--', 'LineWidth', 1.0);

formatPlot( ...
    'Iteration number', ...
    'Residual, |R_i|', ...
    'Newton Residual vs Iteration');

hold off;


%% ================================================================
% LOCAL PLOTTING FUNCTION
% ================================================================

function formatPlot(xText, yText, titleText)

    xlabel(xText, 'FontWeight','bold');
    ylabel(yText, 'FontWeight','bold');

    title(titleText, 'FontWeight','normal');

    grid on;
    box on;

end