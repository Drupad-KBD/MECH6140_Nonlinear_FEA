function tangentStiffness = computeExactTangent(displacement, exactTangentFunction)
% Evaluate the analytical derivative generated from the supplied force law.
tangentStiffness = exactTangentFunction(displacement);
end