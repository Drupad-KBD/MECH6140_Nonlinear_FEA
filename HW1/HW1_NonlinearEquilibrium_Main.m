clear;
clc;

%% USER SETTINGS

% Problem
% hardening spring
% internalForceFunction = @(u) 10*u + 10*u.^3;
% exactTangentFunction  = @(u) 10 + 30*u.^2;
% Softening spring
internalForceFunction = @(u) 12*u - u.^3;
exactTangentFunction  = @(u) 12 - 3*u.^2;

% Solver
finalLoad = 16;
initialDisplacement = 0;
numberOfLoadSteps = 1;
methods = ["exact", "modified", "secant", "numerical"];
residualTolerance = 1e-8;
maximumIterations = 50;
numericalPerturbation = 1e-5;

% Display
displayIterations = false;
showLiveNewtonPlot = true;
livePlotPause = 0.8;

% GIF saving
saveLiveGIF = true;
gifDelay = 0.4;

% Output
outputFile = "Softening_spring_P=12_3_load_step.mat";

%% SETUP

options.finalLoad = finalLoad;
options.initialDisplacement = initialDisplacement;
options.numberOfLoadSteps = numberOfLoadSteps;
options.methods = methods;
options.residualTolerance = residualTolerance;
options.maximumIterations = maximumIterations;
options.numericalPerturbation = numericalPerturbation;

options.displayIterations = displayIterations;
options.showLiveNewtonPlot = showLiveNewtonPlot;
options.livePlotPause = livePlotPause;

options.saveLiveGIF = saveLiveGIF;
options.gifDelay = gifDelay;
options.gifFolder = fullfile(pwd, 'Newton_Animations');

options.outputFile = outputFile;

if saveLiveGIF && ~exist(options.gifFolder, 'dir')
    mkdir(options.gifFolder);
end
%% RUN SOLVER

results = struct();

for methodIndex = 1:numel(methods)

    methodName = methods(methodIndex);

    results.(methodName) = solveScalarNonlinearEquilibrium( ...
        internalForceFunction, ...
        exactTangentFunction, ...
        options, ...
        methodName);

end
%% SAVE RESULTS

save(outputFile, 'results', 'options');

fprintf('\nSaved results to %s\n', outputFile);
%% SUMMARY

fprintf('\nMethod       Status       Final u        Corrections\n');
fprintf('----------------------------------------------------\n');

for methodIndex = 1:numel(methods)

    methodName = methods(methodIndex);
    r = results.(methodName);

    if r.converged
        statusText = 'converged';
    else
        statusText = 'failed';
    end

    fprintf('%-12s %-12s %-14.8g %d\n', ...
        methodName, ...
        statusText, ...
        r.finalDisplacement, ...
        r.numberOfCorrections);
end
