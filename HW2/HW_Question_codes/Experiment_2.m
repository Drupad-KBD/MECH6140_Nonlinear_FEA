%% ================================================================
% MECH 6140 - HW2
% EXPERIMENT 2 - Stretch with Superposed Rotation
%
% Plot W versus rotation angle
% Compare E, S, and sigma at 0 and 45 degrees
% ================================================================

clear;
clc;

format short g;


%% ================================================================
% Inputs
% ================================================================

Nodal_Locations = [ ...
    0 0;
    2 0;
    2 1;
    0 1];

Youngs_Modulus = 70000;
Poissons_Ratio = 0.30;

lambda_L = ...
    (Youngs_Modulus*Poissons_Ratio) / ...
    ((1+Poissons_Ratio)*(1-2*Poissons_Ratio));

mu = ...
    Youngs_Modulus / ...
    (2*(1+Poissons_Ratio));

ShapeFun_derivative = ...
    @Derivative_Of_ShapeFun_4Noded;


%% ================================================================
% Reference coordinates
% ================================================================

x_cap = ...
    X_2_to_X_1converter( ...
    Nodal_Locations);


%% ================================================================
% Select one Gauss point
% ================================================================

[Gauss_location,~] = ...
    gauss_location_weight(2);

xi_cordi = [ ...
    Gauss_location(1), ...
    Gauss_location(1)];


%% ================================================================
% Rotation angles
% ================================================================

theta_values = ...
    0:1:90;

number_theta = ...
    length(theta_values);

W_values = ...
    zeros(number_theta,1);


%% ================================================================
% Rotate stretched element
% ================================================================

for k = 1:number_theta

    theta = ...
        theta_values(k);


    %% ------------------------------------------------------------
    % Stretch first, then rotate
    %
    % x = Q(theta) * [1.2X
    %                  0.9Y]
    % -------------------------------------------------------------

    u = @(X,Y) ...
        cosd(theta)*(1.2*X) ...
        - sind(theta)*(0.9*Y) ...
        - X;

    v = @(X,Y) ...
        sind(theta)*(1.2*X) ...
        + cosd(theta)*(0.9*Y) ...
        - Y;


    %% Prescribed nodal displacement vector

    q = ...
        Prescribed_Nodal_Displacements_cal( ...
        Nodal_Locations, ...
        u,v);


    %% Displacement gradient

    H = ...
        Displacement_Gradient_cal( ...
        ShapeFun_derivative, ...
        xi_cordi, ...
        x_cap, ...
        q);


    %% Deformation gradient

    F = ...
        Deformation_Gradient_cal(H);


    %% Deformation Jacobian

    J = ...
        Deformation_Jacobian_cal(F);


    %% Green-Lagrange strain

    [E_GL,~] = ...
        Green_Lagrange_Strain_cal(H);


    %% Strain energy density

    W = ...
        Strain_Energy_Density_cal( ...
        E_GL, ...
        lambda_L, ...
        mu);


    %% Second Piola-Kirchhoff stress

    S = ...
        Second_Piola_Kirchhoff_Stress_cal( ...
        E_GL, ...
        lambda_L, ...
        mu);


    %% Cauchy stress

    sigma = ...
        Cauchy_Stress_cal( ...
        F,S,J);


    %% Store strain energy

    W_values(k) = ...
        W;


    %% ------------------------------------------------------------
    % Store quantities at 0 and 45 degrees
    % -------------------------------------------------------------

    if theta == 0

        E_0 = E_GL;
        S_0 = S;
        sigma_0 = sigma;

    end


    if theta == 45

        E_45 = E_GL;
        S_45 = S;
        sigma_45 = sigma;

    end

end


%% ================================================================
% Display comparison
% ================================================================

disp(' ');
disp('================================================');
disp('GREEN-LAGRANGE STRAIN');
disp('================================================');

disp('E at 0 degrees');
disp(E_0);

disp('E at 45 degrees');
disp(E_45);


disp(' ');
disp('================================================');
disp('SECOND PIOLA-KIRCHHOFF STRESS');
disp('================================================');

disp('S at 0 degrees');
disp(S_0);

disp('S at 45 degrees');
disp(S_45);


disp(' ');
disp('================================================');
disp('CAUCHY STRESS');
disp('================================================');

disp('sigma at 0 degrees');
disp(sigma_0);

disp('sigma at 45 degrees');
disp(sigma_45);


%% ================================================================
% Strain energy plot
% ================================================================

figure;

plot( ...
    theta_values, ...
    W_values, ...
    'LineWidth',1.5);

xlabel('Rotation Angle, \theta (degrees)');

ylabel('Strain Energy Density, W');

title( ...
    'Stretch with Superposed Rigid Rotation', ...
    'FontWeight','normal');

grid on;
box on;


%% Clean Y-axis

W_mean = mean(W_values);

ylim([ ...
    W_mean - 0.1, ...
    W_mean + 0.1]);

ytickformat('%.2f');


%% ================================================================
% Results table
% ================================================================

Experiment2_Results = table( ...
    theta_values.', ...
    W_values, ...
    'VariableNames',{ ...
    'Rotation_Angle_deg', ...
    'Strain_Energy_Density'});


%% ================================================================
% Export to Excel
% ================================================================

Excel_File = ...
    'HW2_Experiment2_StretchRotation.xlsx';

if isfile(Excel_File)

    delete(Excel_File);

end

writetable( ...
    Experiment2_Results, ...
    Excel_File, ...
    'Sheet', ...
    'Energy vs Rotation');

fprintf('\nResults saved to:\n%s\n',Excel_File);