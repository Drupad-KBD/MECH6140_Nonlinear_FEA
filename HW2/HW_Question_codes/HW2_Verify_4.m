%% ================================================================
% MECH 6140 - HW2
% VERIFY 4 - Resultant Force and Moment
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

thickness = 1;

lambda_L = ...
    (Youngs_Modulus*Poissons_Ratio) / ...
    ((1+Poissons_Ratio)*(1-2*Poissons_Ratio));

mu = ...
    Youngs_Modulus / ...
    (2*(1+Poissons_Ratio));

ShapeFun_derivative = ...
    @Derivative_Of_ShapeFun_4Noded;

gamma = 0.5;


%% ================================================================
% Cases
% ================================================================

cases = [ ...
    "stretch", ...
    "shear"];

Resultant_Force_X = ...
    zeros(length(cases),1);

Resultant_Force_Y = ...
    zeros(length(cases),1);

Resultant_Force_Magnitude = ...
    zeros(length(cases),1);

Resultant_Moment = ...
    zeros(length(cases),1);


%% ================================================================
% Run cases
% ================================================================

for k = 1:length(cases)

    caseName = ...
        cases(k);


    %% ------------------------------------------------------------
    % Prescribed motion
    % -------------------------------------------------------------

    switch caseName

        case "stretch"

            u = @(X,Y)  0.2*X;

            v = @(X,Y) -0.1*Y;


        case "shear"

            u = @(X,Y) gamma*Y;

            v = @(X,Y) 0;

    end


    %% ============================================================
    % Nodal displacement vector
    % ============================================================

    q = ...
        Prescribed_Nodal_Displacements_cal( ...
        Nodal_Locations, ...
        u,v);


    %% ============================================================
    % Current nodal coordinates
    % ============================================================

    Current_Nodal_Locations = ...
        Nodal_Locations;

    Current_Nodal_Locations(:,1) = ...
        Nodal_Locations(:,1) ...
        + q(1:2:end);

    Current_Nodal_Locations(:,2) = ...
        Nodal_Locations(:,2) ...
        + q(2:2:end);


    %% ============================================================
    % Reference coordinates
    % ============================================================

    x_cap = ...
        X_2_to_X_1converter( ...
        Nodal_Locations);


    %% ============================================================
    % Gauss integration
    % ============================================================

    [Gauss_location,Weight] = ...
        gauss_location_weight(2);

    Internal_Force = ...
        zeros(8,1);


    for i = 1:2

        for j = 1:2

            xi_cordi = [ ...
                Gauss_location(i), ...
                Gauss_location(j)];

            Gauss_Weight = ...
                Weight(i)*Weight(j);


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


            %% Green-Lagrange strain

            [E_GL,~] = ...
                Green_Lagrange_Strain_cal(H);


            %% Second Piola-Kirchhoff stress

            S = ...
                Second_Piola_Kirchhoff_Stress_cal( ...
                E_GL, ...
                lambda_L, ...
                mu);


            %% Nonlinear B matrix

            B = ...
                Nonlinear_B_Matrix_cal( ...
                ShapeFun_derivative, ...
                xi_cordi, ...
                x_cap, ...
                F);


            %% Internal force contribution

            f_int_GP = ...
                Internal_Force_GP_cal( ...
                B, ...
                S, ...
                ShapeFun_derivative, ...
                xi_cordi, ...
                x_cap, ...
                Gauss_Weight, ...
                thickness);


            Internal_Force = ...
                Internal_Force ...
                + f_int_GP;

        end

    end


    %% ============================================================
    % Separate x and y nodal forces
    % ============================================================

    f_x = ...
        Internal_Force(1:2:end);

    f_y = ...
        Internal_Force(2:2:end);


    %% ============================================================
    % Resultant internal nodal force
    % ============================================================

    Resultant_Force_X(k) = ...
        sum(f_x);

    Resultant_Force_Y(k) = ...
        sum(f_y);

    Resultant_Force_Magnitude(k) = ...
        sqrt( ...
        Resultant_Force_X(k)^2 ...
        + Resultant_Force_Y(k)^2);


    %% ============================================================
    % Resultant moment
    %
    % Mz = sum(x*Fy - y*Fx)
    % ============================================================

    x = ...
        Current_Nodal_Locations(:,1);

    y = ...
        Current_Nodal_Locations(:,2);

    Resultant_Moment(k) = ...
        sum( ...
        x.*f_y ...
        - y.*f_x);

end


%% ================================================================
% Results
% ================================================================

Case = ...
    cases.';

Verify4_Results = table( ...
    Case, ...
    Resultant_Force_X, ...
    Resultant_Force_Y, ...
    Resultant_Force_Magnitude, ...
    Resultant_Moment);


disp(' ');
disp('================================================');
disp('VERIFY 4 - RESULTANT FORCE AND MOMENT');
disp('================================================');

disp(Verify4_Results);


%% ================================================================
% Export to Excel
% ================================================================

Excel_File = ...
    'HW2_Verify4_ForceMoment.xlsx';

if isfile(Excel_File)

    delete(Excel_File);

end

writetable( ...
    Verify4_Results, ...
    Excel_File, ...
    'Sheet', ...
    'Verify 4');

fprintf('\nResults saved to:\n%s\n',Excel_File);