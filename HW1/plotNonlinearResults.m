function plotNonlinearResults(results, internalForceFunction, plotOptions)
% Plot saved histories without rerunning the nonlinear solver.
methodNames = string(fieldnames(results));
allDisplacements = [];
allLoads = [];
for methodIndex = 1:numel(methodNames)
    methodResult = results.(methodNames(methodIndex));
    allDisplacements = [allDisplacements; methodResult.history.displacement; ...
        methodResult.history.nextDisplacement]; %#ok<AGROW>
    allLoads = [allLoads; methodResult.history.targetLoad]; %#ok<AGROW>
end

finiteDisplacements = allDisplacements(isfinite(allDisplacements));
finiteLoads = allLoads(isfinite(allLoads));
if isempty(finiteDisplacements) || isempty(finiteLoads)
    return;
end

if plotOptions.showForceDisplacement
    figure('Name', 'Force-Displacement and Newton Paths');
    hold on;
    displacementRange = localPlotRange(finiteDisplacements);
    plot(displacementRange, internalForceFunction(displacementRange), ...
        'k-', 'LineWidth', 1.5, 'DisplayName', 'Internal force');
    plot(displacementRange, zeros(size(displacementRange)) + max(finiteLoads), ...
        'k--', 'DisplayName', 'Final target load');

    for methodIndex = 1:numel(methodNames)
        methodName = methodNames(methodIndex);
        methodResult = results.(methodName);
        history = methodResult.history;
        validRows = isfinite(history.displacement) & ...
            isfinite(history.internalForce) & isfinite(history.nextDisplacement);
        if plotOptions.showNewtonIterationPaths
            for rowIndex = find(validRows).'
                nextDisplacement = history.nextDisplacement(rowIndex);
                nextForce = internalForceFunction(nextDisplacement);
                plot([history.displacement(rowIndex), nextDisplacement, nextDisplacement], ...
                    [history.internalForce(rowIndex), history.targetLoad(rowIndex), nextForce], ...
                    '.-', 'DisplayName', char(methodName));
            end
        end
    end
    xlabel('Displacement');
    ylabel('Force');
    title('Nonlinear Equilibrium and Newton Marching Paths');
    grid on;
    legend('Location', 'best');
    hold off;
end

if plotOptions.showResidualHistory
    figure('Name', 'Residual Convergence History');
    hold on;
    for methodIndex = 1:numel(methodNames)
        methodName = methodNames(methodIndex);
        history = results.(methodName).history;
        validRows = isfinite(history.absoluteResidual) & ...
            history.absoluteResidual > 0;
        if plotOptions.showMethodComparison || methodIndex == 1
            semilogy(history.iteration(validRows), ...
                history.absoluteResidual(validRows), '.-', ...
                'DisplayName', char(methodName));
        end
    end
    xlabel('Nonlinear iteration');
    ylabel('|R|');
    title('Residual Convergence');
    grid on;
    legend('Location', 'best');
    hold off;
end

if plotOptions.showCorrectionHistory
    figure('Name', 'Displacement Correction History');
    hold on;
    for methodIndex = 1:numel(methodNames)
        methodName = methodNames(methodIndex);
        history = results.(methodName).history;
        validRows = isfinite(history.absoluteCorrection) & ...
            history.absoluteCorrection > 0;
        if plotOptions.showMethodComparison || methodIndex == 1
            semilogy(history.iteration(validRows), ...
                history.absoluteCorrection(validRows), '.-', ...
                'DisplayName', char(methodName));
        end
    end
    xlabel('Nonlinear iteration');
    ylabel('|delta u|');
    title('Displacement Correction History');
    grid on;
    legend('Location', 'best');
    hold off;
end
end

function plotRange = localPlotRange(values)
minimumValue = min(values);
maximumValue = max(values);
span = maximumValue - minimumValue;
if span == 0
    span = max(1, abs(minimumValue));
end
plotRange = linspace(minimumValue - 0.1 * span, ...
    maximumValue + 0.1 * span, 400);
end