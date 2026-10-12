%% ================================================================
% MECH 6140 - HW2
% EXPERIMENT 3 - Small-Deformation Limit
%
% Compare nonlinear internal force with linear FE prediction Kq
%
% q = alpha*q_bar
% ================================================================

clear;
clc;
close all;

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


%% ================================================================
% Combined stretch + shear displacement field
%
% u_bar = 0.2X + 0.5Y
% v_bar = -0.1Y
% ================================================================

u_bar = @(X,Y) ...
    0.2*X + 0.5*Y;

v_bar = @(X,Y) ...
    -0.1*Y;


%% ================================================================
% Base nodal displacement vector
% ================================================================

q_bar = ...
    Prescribed_Nodal_Displacements_cal( ...
    Nodal_Locations, ...
    u_bar, ...
    v_bar);


%% ================================================================
% Linear FE stiffness matrix
% ================================================================

D_linear = ...
    D_matrix_PlaneStrain( ...
    Youngs_Modulus, ...
    Poissons_Ratio);

K = ...
    Gauss_Element_stiffness_cal( ...
    ShapeFun_derivative, ...
    Nodal_Locations, ...
    D_linear, ...
    2, ...
    thickness);


%% ================================================================
% Reference coordinates
% ================================================================

x_cap = ...
    X_2_to_X_1converter( ...
    Nodal_Locations);


%% ================================================================
% Gauss points
% ================================================================

[Gauss_location,Weight] = ...
    gauss_location_weight(2);


%% ================================================================
% Alpha values
%
% Small deformation --> finite deformation
% ================================================================

alpha_values = ...
    logspace(-4,0.5,25);

number_alpha = ...
    length(alpha_values);


%% ================================================================
% Allocate results
% ================================================================

f_int_Nonlinear_All = ...
    zeros(8,number_alpha);

f_int_Linear_All = ...
    zeros(8,number_alpha);

Nonlinear_Force_Magnitude = ...
    zeros(number_alpha,1);

Linear_Force_Magnitude = ...
    zeros(number_alpha,1);

Relative_Error = ...
    zeros(number_alpha,1);


%% ================================================================
% Loop through alpha values
% ================================================================

for k = 1:number_alpha

    alpha = ...
        alpha_values(k);


    %% ------------------------------------------------------------
    % Scale displacement vector
    %
    % q = alpha*q_bar
    % ------------------------------------------------------------

    q = ...
        alpha*q_bar;


    %% ============================================================
    % Nonlinear internal-force vector
    % ============================================================

    f_int_Nonlinear = ...
        zeros(8,1);


    for i = 1:2

        for j = 1:2

            %% ----------------------------------------------------
            % Gauss point
            % ----------------------------------------------------

            xi_cordi = [ ...
                Gauss_location(i), ...
                Gauss_location(j)];

            Gauss_Weight = ...
                Weight(i)*Weight(j);


            %% ----------------------------------------------------
            % Displacement gradient
            % ----------------------------------------------------

            H = ...
                Displacement_Gradient_cal( ...
                ShapeFun_derivative, ...
                xi_cordi, ...
                x_cap, ...
                q);


            %% ----------------------------------------------------
            % Deformation gradient
            % ----------------------------------------------------

            F = ...
                Deformation_Gradient_cal(H);


            %% ----------------------------------------------------
            % Green-Lagrange strain
            % ----------------------------------------------------

            [E_GL,~] = ...
                Green_Lagrange_Strain_cal(H);


            %% ----------------------------------------------------
            % Second Piola-Kirchhoff stress
            % ----------------------------------------------------

            S = ...
                Second_Piola_Kirchhoff_Stress_cal( ...
                E_GL, ...
                lambda_L, ...
                mu);


            %% ----------------------------------------------------
            % Nonlinear B matrix
            % ----------------------------------------------------

            B = ...
                Nonlinear_B_Matrix_cal( ...
                ShapeFun_derivative, ...
                xi_cordi, ...
                x_cap, ...
                F);


            %% ----------------------------------------------------
            % Gauss-point internal force
            % ----------------------------------------------------

            f_int_GP = ...
                Internal_Force_GP_cal( ...
                B, ...
                S, ...
                ShapeFun_derivative, ...
                xi_cordi, ...
                x_cap, ...
                Gauss_Weight, ...
                thickness);


            %% ----------------------------------------------------
            % Add Gauss-point contribution
            % ----------------------------------------------------

            f_int_Nonlinear = ...
                f_int_Nonlinear ...
                + f_int_GP;

        end

    end


    %% ============================================================
    % Linear FE prediction
    %
    % f_linear = K*q
    % ============================================================

    f_int_Linear = ...
        K*q;


    %% ============================================================
    % Store complete force vectors
    % ============================================================

    f_int_Nonlinear_All(:,k) = ...
        f_int_Nonlinear;

    f_int_Linear_All(:,k) = ...
        f_int_Linear;


    %% ============================================================
    % Force magnitudes
    % ============================================================

    Nonlinear_Force_Magnitude(k) = ...
        norm(f_int_Nonlinear);

    Linear_Force_Magnitude(k) = ...
        norm(f_int_Linear);


    %% ============================================================
    % Relative error
    %
    % ||f_NL - Kq||
    % --------------
    %    ||f_NL||
    % ============================================================

    Relative_Error(k) = ...
        norm( ...
        f_int_Nonlinear ...
        - f_int_Linear) ...
        / ...
        norm(f_int_Nonlinear);

end


%% ================================================================
% Relative error percentage
% ================================================================

Relative_Error_Percent = ...
    100*Relative_Error;


%% ================================================================
% FIGURE 1
%
% Direct comparison of all 8 force-vector components
% ================================================================

figure;

tiledlayout(4,2);


for dof = 1:8

    nexttile;

    semilogx( ...
        alpha_values, ...
        f_int_Nonlinear_All(dof,:), ...
        '-', ...
        'LineWidth',1.5);

    hold on;

    semilogx( ...
        alpha_values, ...
        f_int_Linear_All(dof,:), ...
        '-', ...
        'LineWidth',1.5);


    xlabel('\alpha');

    ylabel( ...
        sprintf('Force DOF %d',dof));

    title( ...
        sprintf('DOF %d',dof), ...
        'FontWeight','normal');

    grid on;
    box on;


    if dof == 1

        legend( ...
            'Nonlinear f_{int}', ...
            'Linear Kq', ...
            'Location','best');

    end

end


sgtitle( ...
    'Nonlinear Internal Force vs Linear Kq', ...
    'FontWeight','normal');


%% ================================================================
% FIGURE 2
%
% Force magnitude comparison
% ================================================================

figure;

loglog( ...
    alpha_values, ...
    Nonlinear_Force_Magnitude, ...
    '-', ...
    'LineWidth',1.5);

hold on;

loglog( ...
    alpha_values, ...
    Linear_Force_Magnitude, ...
    '-', ...
    'LineWidth',1.5);


xlabel('\alpha');

ylabel('Internal Force Magnitude');

title( ...
    'Nonlinear vs Linear Internal Force Magnitude', ...
    'FontWeight','normal');

legend( ...
    '||f_{int}^{NL}||', ...
    '||Kq||', ...
    'Location','northwest');

grid on;
box on;


%% ================================================================
% FIGURE 3
%
% Relative error
%
% ||f_NL - Kq|| / ||f_NL||
% ================================================================

figure;

loglog( ...
    alpha_values, ...
    Relative_Error, ...
    '-', ...
    'LineWidth',1.5);


xlabel('\alpha');

ylabel( ...
    '||f_{int}^{NL} - Kq|| / ||f_{int}^{NL}||');

title( ...
    'Relative Error: Nonlinear vs Linear Internal Force', ...
    'FontWeight','normal');

grid on;
box on;


%% ================================================================
% Results table
% ================================================================

Results = ...
    table( ...
    alpha_values.', ...
    Nonlinear_Force_Magnitude, ...
    Linear_Force_Magnitude, ...
    Relative_Error, ...
    Relative_Error_Percent, ...
    'VariableNames',{ ...
    'Alpha', ...
    'Nonlinear_Force_Magnitude', ...
    'Linear_Force_Magnitude', ...
    'Relative_Error', ...
    'Relative_Error_Percent'});


disp(' ');
disp('================================================');
disp('EXPERIMENT 3 - SMALL-DEFORMATION LIMIT');
disp('================================================');

disp(Results);


%% ================================================================
% Force component tables
% ================================================================

Nonlinear_Table = ...
    array2table( ...
    f_int_Nonlinear_All.', ...
    'VariableNames',{ ...
    'NL_DOF1','NL_DOF2','NL_DOF3','NL_DOF4', ...
    'NL_DOF5','NL_DOF6','NL_DOF7','NL_DOF8'});


Linear_Table = ...
    array2table( ...
    f_int_Linear_All.', ...
    'VariableNames',{ ...
    'Linear_DOF1','Linear_DOF2', ...
    'Linear_DOF3','Linear_DOF4', ...
    'Linear_DOF5','Linear_DOF6', ...
    'Linear_DOF7','Linear_DOF8'});


Force_Comparison = ...
    [ ...
    table(alpha_values.','VariableNames',{'Alpha'}), ...
    Nonlinear_Table, ...
    Linear_Table];


%% ================================================================
% Export to Excel
% ================================================================

Excel_File = ...
    'HW2_Experiment3_SmallDeformation.xlsx';


if isfile(Excel_File)

    delete(Excel_File);

end


writetable( ...
    Results, ...
    Excel_File, ...
    'Sheet', ...
    'Summary');


writetable( ...
    Force_Comparison, ...
    Excel_File, ...
    'Sheet', ...
    'Force Vector Comparison');


fprintf( ...
    '\nResults saved to:\n%s\n', ...
    Excel_File);


figure;

semilogx( ...
    alpha_values, ...
    Relative_Error, ...
    '-', ...
    'LineWidth',1.5);

xlabel('\alpha');

ylabel( ...
    '||f_{int}^{NL} - Kq|| / ||f_{int}^{NL}||');

title( ...
    'Relative Error: Nonlinear vs Linear Internal Force', ...
    'FontWeight','normal');

grid on;
box on;


%% ================================================================
% FIGURE 2
%
% Force magnitude comparison
% ================================================================

figure;

semilogx( ...
    alpha_values, ...
    Nonlinear_Force_Magnitude, ...
    '-', ...
    'LineWidth',1.5);

hold on;

semilogx( ...
    alpha_values, ...
    Linear_Force_Magnitude, ...
    '-', ...
    'LineWidth',1.5);

xlabel('\alpha');

ylabel('Internal Force Magnitude');

title( ...
    'Nonlinear vs Linear Internal Force Magnitude', ...
    'FontWeight','normal');

legend( ...
    '||f_{int}^{NL}||', ...
    '||Kq||', ...
    'Location','northwest');

grid on;
box on;

