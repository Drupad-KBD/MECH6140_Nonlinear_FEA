function diagnostics = checkConvergence( ...
    residual, displacementCorrection, targetLoad, currentDisplacement, ...
    residualTolerance)
diagnostics.absoluteResidual = abs(residual);
if targetLoad == 0
    diagnostics.relativeResidual = NaN;
else
    diagnostics.relativeResidual = diagnostics.absoluteResidual / abs(targetLoad);
end

if isfinite(displacementCorrection)
    diagnostics.absoluteCorrection = abs(displacementCorrection);
    if currentDisplacement == 0
        diagnostics.relativeCorrection = NaN;
    else
        diagnostics.relativeCorrection = diagnostics.absoluteCorrection / ...
            abs(currentDisplacement);
    end
else
    diagnostics.absoluteCorrection = NaN;
    diagnostics.relativeCorrection = NaN;
end

diagnostics.isConverged = diagnostics.absoluteResidual <= residualTolerance;
end