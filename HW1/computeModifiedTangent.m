function tangentStiffness = computeModifiedTangent(displacement, exactTangentFunction)
% Modified Newton uses the tangent evaluated at the load-step starting state.
tangentStiffness = computeExactTangent(displacement, exactTangentFunction);
end