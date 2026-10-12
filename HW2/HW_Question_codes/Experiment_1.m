%% ================================================================
% MECH 6140 - HW2
% EXPERIMENT 1 - Rigid Rotation
%
% Plot ||epsilon|| and ||E|| versus rotation angle
% ================================================================

clear;
clc;

format short g;


%% ================================================================
% Reference nodal coordinates
% ================================================================

Nodal_Locations = [ ...
    0 0;
    2 0;
    2 1;
    0 1];


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


epsilon_norm = ...
    zeros(number_theta,1);

E_GL_norm = ...
    zeros(number_theta,1);


%% ================================================================
% Rotate element
% ================================================================

for k = 1:number_theta

    theta = ...
        theta_values(k);


    %% Rigid rotation displacement field
    %
    % x = cos(theta)X - sin(theta)Y
    % y = sin(theta)X + cos(theta)Y
    %
    % u = x - X
    % v = y - Y

    u = @(X,Y) ...
        cosd(theta)*X ...
        - sind(theta)*Y ...
        - X;

    v = @(X,Y) ...
        sind(theta)*X ...
        + cosd(theta)*Y ...
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


    %% Infinitesimal strain

    [epsilon,~] = ...
        Infinitesimal_Strain_cal(H);


    %% Green-Lagrange strain

    [E_GL,~] = ...
        Green_Lagrange_Strain_cal(H);


    %% Tensor norms

    epsilon_norm(k) = ...
        norm(epsilon,'fro');

    E_GL_norm(k) = ...
        norm(E_GL,'fro');

end


%% ================================================================
% Results table
% ================================================================

Experiment1_Results = table( ...
    theta_values.', ...
    epsilon_norm, ...
    E_GL_norm, ...
    'VariableNames',{ ...
    'Rotation_Angle_deg', ...
    'Infinitesimal_Strain_Norm', ...
    'Green_Lagrange_Strain_Norm'});


disp(' ');
disp('================================================');
disp('EXPERIMENT 1 - RIGID ROTATION');
disp('================================================');

disp(Experiment1_Results);


%% ================================================================
% Plot
% ================================================================

figure;

plot( ...
    theta_values, ...
    epsilon_norm, ...
    'LineWidth',1.5);

hold on;

plot( ...
    theta_values, ...
    E_GL_norm, ...
    'LineWidth',1.5);

xlabel('Rotation Angle, \theta (degrees)');

ylabel('Strain Norm');

title( ...
    'Rigid Rotation: Infinitesimal vs Green-Lagrange Strain', ...
    'FontWeight','normal');

legend( ...
    '||\epsilon||', ...
    '||E||', ...
    'Location','northwest');

grid on;
box on;


%% ================================================================
% Export to Excel
% ================================================================

Excel_File = ...
    'HW2_Experiment1_RigidRotation.xlsx';

if isfile(Excel_File)

    delete(Excel_File);

end

writetable( ...
    Experiment1_Results, ...
    Excel_File, ...
    'Sheet', ...
    'Experiment 1');

fprintf('\nResults saved to:\n%s\n',Excel_File);