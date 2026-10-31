
%% ================================================================
% MECH 6140 - HW3
% VERIFY STRESS TRANSFORMATION AND MATERIAL MODELS
%
% Uses existing OOP classes:
%   fea.materials.LinearElastic
%   fea.materials.LogElastic
%   fea.mechanics.FiniteStrain
%
% No existing class modifications required.
% ================================================================

clearvars;
clc;
close all;

format long g;

%% 1. MATERIAL PROPERTIES

Eyoung = 70000;
nu = 0.30;

StVK = fea.materials.LinearElastic(Eyoung,nu);
LogElastic = fea.materials.LogElastic(Eyoung,nu);

%% 2. DEFORMATION CASES

% Each column gives [lambdaX; lambdaY; shear].
% Simple F matrix:
%
% F = [lambdaX, shear;
%      0,       lambdaY]
%
% Case 1: undeformed
% Case 2: small extension
% Case 3: finite extension
% Case 4: finite extension + transverse contraction
% Case 5: finite deformation + shear

cases = [ ...
    1.00, 1.01, 1.20, 1.30, 1.25;
    1.00, 1.00, 1.00, 0.85, 0.90;
    0.00, 0.00, 0.00, 0.00, 0.20];

caseNames = { ...
    'Reference', ...
    'Small tension', ...
    'Finite tension', ...
    'Tension + contraction', ...
    'Finite deformation + shear'};

nCases = size(cases,2);

%% 3. STORAGE

relativeError = zeros(nCases,1);
materialDifference = zeros(nCases,1);

S11_StVK = zeros(nCases,1);
S11_Recovered = zeros(nCases,1);
S11_LogElastic = zeros(nCases,1);

sigma11_StVK = zeros(nCases,1);
sigma11_LogElastic = zeros(nCases,1);

%% 4. TEST ALL DEFORMATION CASES

fprintf('\n');
fprintf('==========================================================\n');
fprintf(' STVK STRESS TRANSFORMATION VERIFICATION\n');
fprintf('==========================================================\n');

for k = 1:nCases

    %% Deformation gradient

    lambdaX = cases(1,k);
    lambdaY = cases(2,k);
    shear = cases(3,k);

    F = [lambdaX, shear;
         0,       lambdaY];

    J = det(F);

    assert(J>0,'Invalid deformation gradient.');

    H = F-eye(2);

    %% Green-Lagrange strain using existing OOP mechanics

    [EGL,~] = ...
        fea.elements.ElementCalculations. ...
        Green_Lagrange_Strain_cal(H);

    %% A. STVK - DIRECT SECOND PK

    Sdirect = ...
        StVK.Second_Piola_Kirchhoff_Stress_cal(EGL);

    %% B. PUSH FORWARD TO CAUCHY

    sigmaStVK = (F*Sdirect*F.')/J;

    %% C. PULL BACK TO SECOND PK

    Srecovered = ...
        J*(F\sigmaStVK)/F.';

    %% D. LOGELASTIC - SECOND PK FROM EXISTING CLASS

    Slog = ...
        LogElastic.Second_Piola_Kirchhoff_Stress_cal(EGL);

    %% LogElastic Second PK to Cauchy for comparison

    sigmaLog = (F*Slog*F.')/J;

    %% VERIFICATION ERRORS

    relativeError(k) = ...
        norm(Srecovered-Sdirect,'fro') / ...
        max(1,norm(Sdirect,'fro'));

    materialDifference(k) = ...
        norm(Slog-Sdirect,'fro') / ...
        max(1,norm(Sdirect,'fro'));

    S11_StVK(k) = Sdirect(1,1);
    S11_Recovered(k) = Srecovered(1,1);
    S11_LogElastic(k) = Slog(1,1);

    sigma11_StVK(k) = sigmaStVK(1,1);
    sigma11_LogElastic(k) = sigmaLog(1,1);

    %% PRINT RESULTS

    fprintf('\nCase %d: %s\n',k,caseNames{k});
    fprintf('J = %.6f\n',J);

    fprintf('S11 StVK direct       : %.10f\n',Sdirect(1,1));
    fprintf('S11 StVK recovered    : %.10f\n',Srecovered(1,1));
    fprintf('S11 LogElastic        : %.10f\n',Slog(1,1));

    fprintf('Transformation error  : %.6e\n',relativeError(k));
    fprintf('Material difference   : %.6e\n',materialDifference(k));

end

%% 5. SUMMARY TABLE

results = table( ...
    string(caseNames(:)), ...
    S11_StVK, ...
    S11_Recovered, ...
    S11_LogElastic, ...
    relativeError, ...
    materialDifference, ...
    'VariableNames',{ ...
    'Case', ...
    'S11_StVK', ...
    'S11_Recovered', ...
    'S11_LogElastic', ...
    'TransformationError', ...
    'MaterialDifference'});

fprintf('\n========== SUMMARY ==========\n');
disp(results);

%% 6. VERIFICATION

maxError = max(relativeError);

fprintf('\nMaximum stress transformation error: %.6e\n',maxError);

assert(maxError<1e-12, ...
    'StVK stress transformation failed.');

fprintf('PASS: StVK direct and transformed stresses agree.\n');

%% 7. PLOT - STVK TRANSFORMATION

figure('Name','StVK Stress Transformation','Color','w');

plot(1:nCases,S11_StVK,'bo-', ...
    'LineWidth',1.5,'MarkerSize',6);
hold on;

plot(1:nCases,S11_Recovered,'rx--', ...
    'LineWidth',1.2,'MarkerSize',7);

xticks(1:nCases);
xticklabels(caseNames);
xtickangle(20);

xlabel('Deformation case');
ylabel('Second PK Stress S_{11}');

legend('StVK Direct','StVK Transformed', ...
    'Location','best');

title('Stress Transformation Verification', ...
    'FontWeight','normal');

grid on;
box on;

%% 8. PLOT - CONSTITUTIVE MATERIAL COMPARISON

figure('Name','StVK vs LogElastic','Color','w');

plot(1:nCases,S11_StVK,'bo-', ...
    'LineWidth',1.5,'MarkerSize',6);
hold on;

plot(1:nCases,S11_LogElastic,'rs--', ...
    'LineWidth',1.3,'MarkerSize',6);

xticks(1:nCases);
xticklabels(caseNames);
xtickangle(20);

xlabel('Deformation case');
ylabel('Second PK Stress S_{11}');

legend('StVK','LogElastic', ...
    'Location','best');

title('Different Constitutive Laws', ...
    'FontWeight','normal');

grid on;
box on;
