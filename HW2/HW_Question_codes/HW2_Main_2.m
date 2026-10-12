%% ================================================================
% MECH 6140 - HW2
% Run deformation cases and export to Excel
% ================================================================

clear; clc;

format short g;


%% Inputs

Nodal_Locations = [ ...
    0 0;
    2 0;
    2 1;
    0 1];

E  = 70000;
nu = 0.3;

lambda_L = ...
    (E*nu)/((1+nu)*(1-2*nu));

mu = ...
    E/(2*(1+nu));

thickness = 1;

ShapeFun_derivative = ...
    @Derivative_Of_ShapeFun_4Noded;

theta = 30;      % rotation angle

gamma = 0.5;     % simple shear


%% ---------------------------------------------------------------
% SELECT CASE
%
% "all"
% "translation"
% "rotation"
% "stretch"
% "shear"
% "stretch_rotation"
% ---------------------------------------------------------------

caseToRun = "all";

allCases = [ ...
    "translation", ...
    "rotation", ...
    "stretch", ...
    "shear", ...
    "stretch_rotation"];

if caseToRun == "all"

    cases = allCases;

else

    cases = caseToRun;

end


%% Excel file

Excel_File = ...
    'HW2_data.xlsx';

if isfile(Excel_File)

    delete(Excel_File);

end


%% ================================================================
% Internal force storage for all cases
% ================================================================

Internal_Force_All = ...
    zeros(length(cases),8);


%% ================================================================
% Run selected cases
% ================================================================

for k = 1:length(cases)

    caseName = ...
        cases(k);


    %% Prescribed motion

    switch caseName


        case "translation"

            u = @(X,Y) 0.25;

            v = @(X,Y) -0.10;


        case "rotation"

            u = @(X,Y) ...
                cosd(theta)*X ...
                - sind(theta)*Y ...
                - X;

            v = @(X,Y) ...
                sind(theta)*X ...
                + cosd(theta)*Y ...
                - Y;


        case "stretch"

            u = @(X,Y) 0.2*X;

            v = @(X,Y) -0.1*Y;


        case "shear"

            u = @(X,Y) gamma*Y;

            v = @(X,Y) 0;


        case "stretch_rotation"

            u = @(X,Y) ...
                cosd(theta)*(1.2*X) ...
                - sind(theta)*(0.9*Y) ...
                - X;

            v = @(X,Y) ...
                sind(theta)*(1.2*X) ...
                + cosd(theta)*(0.9*Y) ...
                - Y;

    end


    %% Run calculation

    [Results,Internal_Force] = ...
        HW2_Case_Calculator( ...
        Nodal_Locations, ...
        u,v, ...
        ShapeFun_derivative, ...
        lambda_L,mu, ...
        thickness);


    %% Save Gauss point results to Excel

    writetable( ...
        Results, ...
        Excel_File, ...
        'Sheet', ...
        char(caseName));


    %% Store total internal force

    Internal_Force_All(k,:) = ...
        Internal_Force.';


    %% Display

    fprintf('\n================================================\n');
    fprintf('%s\n',upper(char(caseName)));
    fprintf('================================================\n');

    fprintf('Internal Force Vector:\n');

    disp(Internal_Force);


    fprintf('Finished: %s\n',caseName);

end


%% ================================================================
% Internal force summary
% ================================================================

Case = ...
    cases.';

Internal_Force_Results = table( ...
    Case, ...
    Internal_Force_All(:,1), ...
    Internal_Force_All(:,2), ...
    Internal_Force_All(:,3), ...
    Internal_Force_All(:,4), ...
    Internal_Force_All(:,5), ...
    Internal_Force_All(:,6), ...
    Internal_Force_All(:,7), ...
    Internal_Force_All(:,8), ...
    'VariableNames',{ ...
    'Case', ...
    'f1x','f1y', ...
    'f2x','f2y', ...
    'f3x','f3y', ...
    'f4x','f4y'});


%% Save internal forces to Excel

writetable( ...
    Internal_Force_Results, ...
    Excel_File, ...
    'Sheet', ...
    'Internal_Forces');


%% Display summary

disp(' ');
disp('================================================');
disp('INTERNAL FORCE SUMMARY');
disp('================================================');

disp(Internal_Force_Results);


fprintf('\nResults saved to:\n%s\n',Excel_File);



%% ================================================================
% CASE CALCULATOR
% ================================================================

function [Results,Internal_Force] = ...
    HW2_Case_Calculator( ...
    Nodal_Locations,u,v, ...
    ShapeFun_derivative, ...
    lambda_L,mu, ...
    thickness)


%% Prescribed nodal displacement vector

q = ...
    Prescribed_Nodal_Displacements_cal( ...
    Nodal_Locations,u,v);


%% Reference nodal coordinates

x_cap = ...
    X_2_to_X_1converter( ...
    Nodal_Locations);


%% Gauss points

[gp,w] = ...
    gauss_location_weight(2);


%% Initialize total internal force

Internal_Force = ...
    zeros(8,1);


counter = 0;


%% ================================================================
% Gauss point loop
% ================================================================

for i = 1:2

    for j = 1:2

        counter = ...
            counter + 1;

        xi = [ ...
            gp(i), ...
            gp(j)];

        weight = ...
            w(i)*w(j);


        %% Displacement gradient

        H = ...
            Displacement_Gradient_cal( ...
            ShapeFun_derivative, ...
            xi, ...
            x_cap, ...
            q);


        %% Deformation gradient

        F = ...
            Deformation_Gradient_cal(H);


        %% Deformation Jacobian

        J = ...
            Deformation_Jacobian_cal(F);


        %% Infinitesimal strain

        [eps,~] = ...
            Infinitesimal_Strain_cal(H);


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


        %% --------------------------------------------------------
        % Nonlinear B matrix
        % ---------------------------------------------------------

        B = ...
            Nonlinear_B_Matrix_cal( ...
            ShapeFun_derivative, ...
            xi, ...
            x_cap, ...
            F);


        %% --------------------------------------------------------
        % Internal force contribution at Gauss point
        % ---------------------------------------------------------

        f_int_GP = ...
            Internal_Force_GP_cal( ...
            B, ...
            S, ...
            ShapeFun_derivative, ...
            xi, ...
            x_cap, ...
            weight, ...
            thickness);


        %% --------------------------------------------------------
        % Add Gauss point contribution
        % ---------------------------------------------------------

        Internal_Force = ...
            Internal_Force ...
            + f_int_GP;


        %% ========================================================
        % Store original Gauss point results
        % ========================================================

        GP(counter,1) = ...
            counter;

        Xi(counter,1) = ...
            xi(1);

        Eta(counter,1) = ...
            xi(2);


        H11(counter,1) = H(1,1);

        H12(counter,1) = H(1,2);

        H21(counter,1) = H(2,1);

        H22(counter,1) = H(2,2);


        F11(counter,1) = F(1,1);

        F12(counter,1) = F(1,2);

        F21(counter,1) = F(2,1);

        F22(counter,1) = F(2,2);


        Jdet(counter,1) = ...
            J;


        eps11(counter,1) = eps(1,1);

        eps22(counter,1) = eps(2,2);

        eps12(counter,1) = eps(1,2);


        E11(counter,1) = E_GL(1,1);

        E22(counter,1) = E_GL(2,2);

        E12(counter,1) = E_GL(1,2);


        Energy(counter,1) = ...
            W;


        S11(counter,1) = S(1,1);

        S22(counter,1) = S(2,2);

        S12(counter,1) = S(1,2);


        sig11(counter,1) = sigma(1,1);

        sig22(counter,1) = sigma(2,2);

        sig12(counter,1) = sigma(1,2);

    end

end


%% ================================================================
% Original Gauss point results table
% ================================================================

Results = table( ...
    GP,Xi,Eta, ...
    H11,H12,H21,H22, ...
    F11,F12,F21,F22, ...
    Jdet, ...
    eps11,eps22,eps12, ...
    E11,E22,E12, ...
    Energy, ...
    S11,S22,S12, ...
    sig11,sig22,sig12);


disp(Results);

end


