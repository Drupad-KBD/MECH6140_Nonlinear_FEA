function tangentStiffness = computeNumericalTangent( ...
    displacement, internalForceFunction, numericalPerturbation)
% Centered finite difference of the internal-force function.
positiveForce = internalForceFunction(displacement + numericalPerturbation);
negativeForce = internalForceFunction(displacement - numericalPerturbation);
tangentStiffness = (positiveForce - negativeForce) / ...
    (2 * numericalPerturbation);
end