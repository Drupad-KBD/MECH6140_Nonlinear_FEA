function livePlotState = updateLiveNewtonPlot( ...
    internalForceFunction, currentDisplacement, currentInternalForce, ...
    targetLoad, displacementCorrection, tangentStiffness, methodName, ...
    loadStep, iteration, livePlotState, options)
% Update one figure with the current Newton construction.
if isempty(livePlotState) || ~isfield(livePlotState, 'figure') || ...
        ~isgraphics(livePlotState.figure)
    livePlotState.figure = figure('Name', 'Live Newton Iteration');
    livePlotState.previousDisplacements = [];
    livePlotState.previousForces = [];
else
    figure(livePlotState.figure);
end

livePlotState.previousDisplacements(end + 1) = currentDisplacement;
livePlotState.previousForces(end + 1) = currentInternalForce;
nextDisplacement = currentDisplacement + displacementCorrection;
plotDisplacements = [livePlotState.previousDisplacements, nextDisplacement];
plotDisplacements = plotDisplacements(isfinite(plotDisplacements));
plotRange = localPlotRange(plotDisplacements);

cla(livePlotState.figure);
hold on;
plot(plotRange, internalForceFunction(plotRange), 'k-', ...
    'LineWidth', 1.5, 'DisplayName', 'Internal force');
plot(plotRange, zeros(size(plotRange)) + targetLoad, 'k--', ...
    'DisplayName', 'Target load');
plot(livePlotState.previousDisplacements, livePlotState.previousForces, ...
    'ko', 'MarkerFaceColor', 'w', 'DisplayName', 'Previous points');
plot(currentDisplacement, currentInternalForce, 'bo', ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Current point');
plot([currentDisplacement, nextDisplacement], ...
    [currentInternalForce, targetLoad], 'r-', ...
    'LineWidth', 1.2, 'DisplayName', 'Tangent step');
plot(nextDisplacement, targetLoad, 'rx', 'LineWidth', 1.5, ...
    'DisplayName', 'Predicted intersection');
plot([nextDisplacement, nextDisplacement], ...
    [targetLoad, internalForceFunction(nextDisplacement)], 'g-', ...
    'DisplayName', 'Return to curve');
xlabel('Displacement');
ylabel('Force');
title(sprintf('%s Newton: Load Step %d, Iteration %d', ...
    methodName, loadStep, iteration));
grid on;
legend('Location', 'best');
hold off;
drawnow;

if options.saveLiveGIF
    frame = getframe(livePlotState.figure);
    [indexedImage, colorMap] = rgb2ind(frame2im(frame), 256);
    gifName = sprintf('Newton_%s_%d_LoadSteps.gif', ...
        char(methodName), options.numberOfLoadSteps);
    gifFile = fullfile(options.gifFolder, gifName);
    if loadStep == 1 && iteration == 1
        imwrite(indexedImage, colorMap, gifFile, 'gif', ...
            'LoopCount', inf, 'DelayTime', options.gifDelay);
    else
        imwrite(indexedImage, colorMap, gifFile, 'gif', ...
            'WriteMode', 'append', 'DelayTime', options.gifDelay);
    end
end

if options.livePlotPause > 0
    pause(options.livePlotPause);
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