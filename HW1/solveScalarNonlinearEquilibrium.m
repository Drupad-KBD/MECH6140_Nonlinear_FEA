function methodResult = solveScalarNonlinearEquilibrium( ...
    internalForceFunction, exactTangentFunction, options, methodName)
% Solve every selected tangent method with the same load-step/Newton algorithm.

numberOfLoadSteps = options.numberOfLoadSteps;
targetLoads = linspace(options.finalLoad / numberOfLoadSteps, ...
    options.finalLoad, numberOfLoadSteps);
convergedDisplacement = options.initialDisplacement;
livePlotState = [];

history = struct('loadStep', [], 'targetLoad', [], 'iteration', [], ...
    'displacement', [], 'internalForce', [], 'residual', [], ...
    'tangent', [], 'displacementCorrection', [], 'nextDisplacement', [], ...
    'absoluteResidual', [], 'relativeResidual', [], ...
    'absoluteCorrection', [], 'relativeCorrection', [], ...
    'isConverged', []);
loadStepSummary = struct('loadStep', [], 'targetLoad', [], ...
    'startingDisplacement', [], 'finalDisplacement', [], ...
    'converged', [], 'iterations', [], 'corrections', []);

historyRow = 0;
summaryRow = 0;
methodFailed = false;
failedLoadStep = NaN;

for loadStep = 1:numberOfLoadSteps
    targetLoad = targetLoads(loadStep);
    startingDisplacement = convergedDisplacement;
    currentDisplacement = startingDisplacement;
    previousDisplacement = NaN;
    previousInternalForce = NaN;
    displacementCorrection = NaN;

    initialTangentStiffness = computeExactTangent( ...
        startingDisplacement, exactTangentFunction);
    if methodName == "modified"
        fixedTangentStiffness = computeModifiedTangent( ...
            startingDisplacement, exactTangentFunction);
    else
        fixedTangentStiffness = NaN;
    end

    if options.displayIterations
        fprintf('\nMethod: %s, Load Step: %d, Target Load: %.10g\n', ...
            methodName, loadStep, targetLoad);
        fprintf(['Load Step | Iter | u | fint | Residual | Tangent | ', ...
            'deltaU | Relative Residual | Relative Correction\n']);
        fprintf(['--------------------------------------------------------------------------------------------------\n']);
    end

    stepConverged = false;
    numberOfCorrections = 0;
    for iteration = 1:options.maximumIterations
        internalForce = internalForceFunction(currentDisplacement);
        residual = targetLoad - internalForce;
        diagnostics = checkConvergence( ...
            residual, NaN, targetLoad, currentDisplacement, ...
            options.residualTolerance);
        invalidIteration = ~isfinite(currentDisplacement) || ...
            ~isfinite(internalForce) || ~isfinite(residual);

        tangentStiffness = NaN;
        nextDisplacement = currentDisplacement;
        displacementCorrection = NaN;
        if ~diagnostics.isConverged && ~invalidIteration
            switch methodName
                case "exact"
                    tangentStiffness = computeExactTangent( ...
                        currentDisplacement, exactTangentFunction);
                case "modified"
                    tangentStiffness = fixedTangentStiffness;
                case "secant"
                    if iteration == 1
                        tangentStiffness = initialTangentStiffness;
                    else
                        tangentStiffness = computeSecantTangent( ...
                            currentDisplacement, previousDisplacement, ...
                            internalForce, previousInternalForce, ...
                            initialTangentStiffness);
                    end
                case "numerical"
                    tangentStiffness = computeNumericalTangent( ...
                        currentDisplacement, internalForceFunction, ...
                        options.numericalPerturbation);
                otherwise
                    error('solveScalarNonlinearEquilibrium:UnknownMethod', ...
                        'Unknown tangent method: %s', methodName);
            end

            if ~isfinite(tangentStiffness) || abs(tangentStiffness) <= eps
                invalidIteration = true;
            else
                displacementCorrection = residual / tangentStiffness;
                nextDisplacement = currentDisplacement + displacementCorrection;
                numberOfCorrections = numberOfCorrections + 1;
                invalidIteration = ~isfinite(displacementCorrection) || ...
                    ~isfinite(nextDisplacement);
            end
            diagnostics = checkConvergence(residual, displacementCorrection, ...
                targetLoad, currentDisplacement, options.residualTolerance);
        end

        historyRow = historyRow + 1;
        history.loadStep(historyRow, 1) = loadStep;
        history.targetLoad(historyRow, 1) = targetLoad;
        history.iteration(historyRow, 1) = iteration;
        history.displacement(historyRow, 1) = currentDisplacement;
        history.internalForce(historyRow, 1) = internalForce;
        history.residual(historyRow, 1) = residual;
        history.tangent(historyRow, 1) = tangentStiffness;
        history.displacementCorrection(historyRow, 1) = displacementCorrection;
        history.nextDisplacement(historyRow, 1) = nextDisplacement;
        history.absoluteResidual(historyRow, 1) = diagnostics.absoluteResidual;
        history.relativeResidual(historyRow, 1) = diagnostics.relativeResidual;
        history.absoluteCorrection(historyRow, 1) = diagnostics.absoluteCorrection;
        history.relativeCorrection(historyRow, 1) = diagnostics.relativeCorrection;
        history.isConverged(historyRow, 1) = diagnostics.isConverged;

        if options.displayIterations
            fprintf('%9d | %4d | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g\n', ...
                loadStep, iteration, currentDisplacement, internalForce, ...
                residual, tangentStiffness, displacementCorrection, ...
                diagnostics.relativeResidual, diagnostics.relativeCorrection);
        end

            if options.showLiveNewtonPlot && ~diagnostics.isConverged && ...
                ~invalidIteration
                livePlotState = updateLiveNewtonPlot( ...
                internalForceFunction, currentDisplacement, internalForce, ...
                targetLoad, displacementCorrection, tangentStiffness, ...
                methodName, loadStep, iteration, livePlotState, options);
            end

        if diagnostics.isConverged
            stepConverged = true;
            convergedDisplacement = currentDisplacement;
            break;
        end

        if invalidIteration
            break;
        end

        previousDisplacement = currentDisplacement;
        previousInternalForce = internalForce;
        currentDisplacement = nextDisplacement;
    end

    summaryRow = summaryRow + 1;
    loadStepSummary.loadStep(summaryRow, 1) = loadStep;
    loadStepSummary.targetLoad(summaryRow, 1) = targetLoad;
    loadStepSummary.startingDisplacement(summaryRow, 1) = startingDisplacement;
    loadStepSummary.finalDisplacement(summaryRow, 1) = convergedDisplacement;
    loadStepSummary.converged(summaryRow, 1) = stepConverged;
    loadStepSummary.iterations(summaryRow, 1) = iteration;
    loadStepSummary.corrections(summaryRow, 1) = numberOfCorrections;

    if stepConverged
        fprintf('Load step %d converged after %d corrections.\n', ...
            loadStep, numberOfCorrections);
    else
        fprintf('Load step %d FAILED after %d iterations; the next load step was not attempted.\n', ...
            loadStep, iteration);
        methodFailed = true;
        failedLoadStep = loadStep;
        break;
    end
end

if isempty(history.iteration)
    finalAbsoluteResidual = NaN;
    finalDisplacement = convergedDisplacement;
    numberOfIterations = 0;
else
    lastHistoryRow = numel(history.iteration);
    finalAbsoluteResidual = history.absoluteResidual(lastHistoryRow);
    finalDisplacement = history.displacement(lastHistoryRow);
    numberOfIterations = lastHistoryRow;
end

methodResult.methodName = methodName;
methodResult.converged = ~methodFailed && all(loadStepSummary.converged);
methodResult.failed = methodFailed;
methodResult.failedLoadStep = failedLoadStep;
methodResult.finalDisplacement = finalDisplacement;
methodResult.finalLoad = targetLoads(min(numel(loadStepSummary.targetLoad), numberOfLoadSteps));
methodResult.numberOfIterations = numberOfIterations;
methodResult.numberOfCorrections = sum(loadStepSummary.corrections);
methodResult.finalAbsoluteResidual = finalAbsoluteResidual;
methodResult.history = history;
methodResult.loadStepSummary = loadStepSummary;
%methodResult.livePlotState = livePlotState;
end