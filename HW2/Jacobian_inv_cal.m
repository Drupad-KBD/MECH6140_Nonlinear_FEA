function [J_inv] = Jacobian_inv_cal(ShapeFun_selection,xi_cordi_4N, x_coordinates)
%% jacobian inverse diagonal matrix

%Calling derivative of shape function
J = Jacobian_cal(ShapeFun_selection,xi_cordi_4N, x_coordinates);

% computing the J inverse
J_Inverse = inv(J);

% storing the J inverse in dioganal format %
J_inv = [J_Inverse zeros(2); zeros(2) J_Inverse];




end